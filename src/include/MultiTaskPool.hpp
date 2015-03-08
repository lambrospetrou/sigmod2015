#ifndef __LP_MULTITASKPOOL__
#define __LP_MULTITASKPOOL__

#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <functional>
#include <queue>
#include <utility>
//#include <iostream>

class MultiTaskPool {
    private:

        const static uint64_t TASKS_FREE_THRESHOLD = 1<<14;

        //template<typename T> using LPFunc = std::function<void(uint32_t, uint32_t, T)>;
        typedef std::function<void(uint32_t, uint32_t, void*)> LPFunc;

        struct Task {
            LPFunc func;
            void* args;
        };

    public:
        MultiTaskPool(uint32_t tsz) : mNumOfThreads(tsz), 
        mWaiting(0), mPoolStopped(false), mPoolRunning(false), mHelperId(tsz), mTasksDeqIdx(0) {
            mTasks.reserve(1<<20);
        }

        void initThreads() {
            for (uint32_t i=0; i<mNumOfThreads; ++i) {
                mThreads.push_back(std::thread(&MultiTaskPool::worker, this, i));
            }
        }

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
            lk.unlock();
            mCondActive.notify_all();
            mCondMaster.notify_all();
            for(auto& t : mThreads) t.join();
        }

        inline void addTask(LPFunc f, void* args) {
            Task t;
            t.func = f;
            t.args = args;
            std::lock_guard<std::mutex> lk(mMutex);
            //mTasks.push(std::move(t));
            mTasks.push_back(std::move(t));
            mCondActive.notify_one();
        }
        inline void addTaskUnsafe(LPFunc f, void* args) {
            Task t;
            t.func = f;
            t.args = args;
            //mTasks.push(std::move(t));
            mTasks.push_back(std::move(t));
        }

    private:
        uint64_t mNumOfThreads;
        std::vector<std::thread> mThreads;

        std::condition_variable mCondActive;
        std::condition_variable mCondMaster;
        std::mutex mMutex;
        uint64_t mWaiting;
        bool mPoolStopped, mPoolRunning;

        std::atomic<uint32_t> mHelperId;
        
        //std::queue<Task> mTasks;
        std::vector<Task> mTasks;
        uint64_t mTasksDeqIdx;


        inline bool isTasksEmpty() { /*return mTasks.empty();*/ return mTasks.size() == mTasksDeqIdx; }

        bool inline canProceed() {
            return ((mPoolRunning && !isTasksEmpty()) || mPoolStopped);
        }

        void worker(uint32_t tid) {
            //cerr << tid << endl;
            while(true) {
                std::unique_lock<std::mutex> lk(mMutex);
                ++mWaiting;
                //cerr << ">" << endl;
                // signal for synchronization!!! - all threads are waiting
                if (mWaiting == mNumOfThreads && isTasksEmpty()) mCondMaster.notify_one();
                if (!canProceed()) mCondActive.wait(lk, [this]{return canProceed();});
                --mWaiting;
                if (mPoolStopped) { lk.unlock(); return; }

                //auto cTask = mTasks.front();
                //mTasks.pop();
                auto cTask = mTasks[mTasksDeqIdx++];
                lk.unlock();
                //cerr << "-" << endl;

                // run the function passing the thread ID
                //cTask.func(mNumOfThreads, tid, std::move(cTask.args));
                cTask.func(mNumOfThreads, tid, cTask.args);
                // wait after the execution for the other threads
                //cerr << "2" << endl;
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
