#ifndef __LP_MULTISINGLETASKPOOL__
#define __LP_MULTISINGLETASKPOOL__

#include "LPThreadpool.hpp"

#include <vector>
#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <functional>
#include <utility>
//#include <iostream>

// This is the Threadpol that combines the functionality of the SingleTaskPool and the MultiTaskPool

class MSTaskPool : public ISingleTaskPool, public IMultiTaskPool {
    private:

        const static uint64_t TASKS_FREE_THRESHOLD = 1<<14;

        //template<typename T> using LPFunc = std::function<void(uint32_t, uint32_t, T)>;
        typedef std::function<void(uint32_t, uint32_t, void*)> LPFunc;

        struct Task {
            LPFunc func;
            void* args;
            Task() : args(nullptr) {}
        };

    public:
        MSTaskPool(uint32_t tsz) : ISingleTaskPool(), IMultiTaskPool(), mNumOfThreads(tsz), 
        mWaiting(0), mPoolStopped(false), mPoolRunning(false), mHelperId(tsz), mTasksDeqIdx(0) {
            mTasks.reserve(1<<16);
            // Single Task funcitonality
            mMaskAll = (static_cast<uint64_t>(1) << mNumOfThreads) - 1;
            mThreadsActivity = 0;
            mSingleExists = false;
        }

        void initThreads() {
            for (uint32_t i=0; i<mNumOfThreads; ++i) {
                mThreads.push_back(std::thread(&MultiTaskPool::worker, this, i));
            }
        }

        // should be matched with a waitAll() in order to ensure that the tasks queue is empty
        inline void startAll() {
            std::unique_lock<std::mutex> lk(mMutex);
            mPoolRunning = true;
            mCondActive.notify_all();
            lk.unlock();
        }

        inline void waitAll() {
            std::unique_lock<std::mutex> lk(mMutex);
            mCondMaster.wait(lk, [this]{
                    return (mWaiting == mNumOfThreads && isTasksEmpty());
                    });
            if (mTasksDeqIdx >= TASKS_FREE_THRESHOLD) {
                //std::cerr << "emptying " << mTasks.size() << ":" << mTasksDeqIdx << std::endl;
                mTasks.resize(0); mTasksDeqIdx = 0;
            }
            lk.unlock();
        }

        inline void helpExecution() {
            uint32_t tid = mHelperId++;
            workerOutside(tid);
            --mHelperId;
        }

        inline void destroy() {
            std::unique_lock<std::mutex> lk(mMutex);
            mPoolStopped = true;
            mCondActive.notify_all();
            mCondMaster.notify_all();
            lk.unlock();
            for(auto& t : mThreads) t.join();
        }

        inline void addTask(LPFunc f, void* args) {
            Task t;
            t.func = f;
            t.args = args;
            std::lock_guard<std::mutex> lk(mMutex);
            mTasks.push_back(std::move(t));
            mCondActive.notify_one();
        }
        inline void addTaskUnsafe(LPFunc f, void* args) {
            Task t;
            t.func = f;
            t.args = args;
            mTasks.push_back(std::move(t));
        }

        ///////////////////////////////////////////////
        // SingleTask functionality
        ///////////////////////////////////////////////

        // THESE TWO CALLS SHOULD BE MATCHED AS A PAIR EVERY TIME THEY ARE CALLED
        // WHENEVER THERE IS A startSingleAll() YOU SHOULD 'NOT' USE THE OTHER startAll()/waitAll()
        // UNTIL YOU CALLED waitSingleAll();

        void startSingleAll(LPFunc f, void *args = nullptr) {
            std::unique_lock<std::mutex> lk(mMutex);
            mSingleTask.func = f;
            mSingleTask.args = args;
            mSingleExists = true;
            //mThreadsActivity = 0; // initializes this in order to identify all the threads 
            mCondActive.notify_all();
            lk.unlock();
        }
        void waitSingleAll() {
            std::unique_lock<std::mutex> lk(mMutex);
            //std::cerr << "mw" << std::endl;
            mCondMaster.wait(lk, [this]{
                    return (mThreadsActivity == mMaskAll) && (mWaiting == mNumOfThreads);
                    });
            // threads are waiting for this to become zero after they execute the function
            //std::cerr << "ml" << std::endl;
            mSingleExists = false;
            mThreadsActivity = 0;
            lk.unlock();
        }

    private:
        uint64_t mNumOfThreads;
        std::vector<std::thread> mThreads;

        std::condition_variable mCondActive;
        std::condition_variable mCondMaster;
        std::mutex mMutex;
        //uint64_t mWaiting;
        std::atomic<uint64_t> mWaiting;
        bool mPoolStopped, mPoolRunning;

        std::atomic<uint32_t> mHelperId;
        
        std::vector<Task> mTasks;
        uint64_t mTasksDeqIdx;

        // required for the SingleTask functionality
        uint64_t mMaskAll;
        uint64_t mThreadsActivity;
        bool mSingleExists;
        Task mSingleTask;

        inline bool isTasksEmpty() { /*return mTasks.empty();*/ return mTasks.size() == mTasksDeqIdx; }

        bool inline canProceed() {
            return ((mPoolRunning && !isTasksEmpty()) || mPoolStopped);
        }
        bool inline canProceed(uint64_t tMask) {
            return (canProceed() || (mSingleExists && (tMask&mThreadsActivity) == 0));
        }

        void worker(uint32_t tid) {
            //cerr << tid << endl;
            uint64_t tMask = static_cast<uint64_t>(1)<<tid;
            while(true) {
                ++mWaiting;
                std::unique_lock<std::mutex> lk(mMutex);
                // signal for synchronization!!! - all threads are waiting
                //std::cerr << "bw" << (tid+2);
                //++mWaiting;
                if (mWaiting == mNumOfThreads && isTasksEmpty()) mCondMaster.notify_one();
                if (!canProceed(tMask)) mCondActive.wait(lk, [tMask, this]{return canProceed(tMask);});
                //--mWaiting;
                //std::cerr << "aw" << (tid+2);
                if (mPoolStopped) { lk.unlock(); return; }

                // priority to the single task
                if (mSingleExists) {
                    mThreadsActivity |= tMask; // required by the SingleTask
                    lk.unlock();
                    --mWaiting;
                    mSingleTask.func(mNumOfThreads, tid, mSingleTask.args);
                } else {
                    auto cTask = mTasks[mTasksDeqIdx++];
                    lk.unlock();
                    --mWaiting;
                    cTask.func(mNumOfThreads, tid, cTask.args); 
                }
            }
        }

        void workerOutside (uint32_t tid) {
            //std::cerr << tid << std::endl;
            while(true) {
                Task cTask; 
                {
                    std::lock_guard<std::mutex> lk(mMutex);
                    //std::unique_lock<std::mutex> lk(mMutex);
                    //std::cerr << ">" << std::endl;
                    // all threads are waiting and no jobs are pending
                    if (mWaiting == mNumOfThreads ||  isTasksEmpty()) {
                        //lk.unlock();
                        return;
                    }
                    if (mPoolStopped || !mPoolRunning) { /*lk.unlock();*/ return; }

                    //cTask = mTasks.front();
                    //mTasks.pop();
                    cTask = mTasks[mTasksDeqIdx++];
                }// end of lock_guard

                //std::cerr << "got job!" << std::endl;
                // run the function passing the thread ID
                cTask.func(mNumOfThreads, tid, cTask.args);
            }
        }

};

#endif
