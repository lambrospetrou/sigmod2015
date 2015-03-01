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

        //template<typename T> using LPFunc = std::function<void(uint32_t, uint32_t, T)>;
        typedef std::function<void(uint32_t, uint32_t, void*)> LPFunc;

        struct Task {
            LPFunc func;
            void* args;
        };

    public:
        MultiTaskPool(uint32_t tsz) : mNumOfThreads(tsz), 
            mWaiting(0), mPoolStopped(false), mPoolRunning(false), mHelperId(tsz) {
        }

        void initThreads() {
            for (uint32_t i=0; i<mNumOfThreads; ++i) {
                mThreads.push_back(std::thread(&MultiTaskPool::worker, this, i));
            }
        }

        void startAll() {
            std::unique_lock<std::mutex> lk(mMutex);
            mPoolRunning = true;
            lk.unlock();
            mCondActive.notify_all();
        }

        void waitAllAndStop() {
            std::unique_lock<std::mutex> lk(mMutex);
            mCondMaster.wait(lk, [this]{
                    return (mWaiting == mNumOfThreads && mTasks.empty());
                    });
            mPoolRunning = false;
            lk.unlock();
        }
        void waitAll() {
            std::unique_lock<std::mutex> lk(mMutex);
            mCondMaster.wait(lk, [this]{
                    return (mWaiting == mNumOfThreads && mTasks.empty());
                    });
            lk.unlock();
        }

        void helpExecution() {
            uint32_t tid = mHelperId++;
            workerOutside(tid);
            --mHelperId;
        }

        void destroy() {
            std::unique_lock<std::mutex> lk(mMutex);
            mPoolStopped = true;
            lk.unlock();
            mCondActive.notify_all();
            mCondMaster.notify_all();
            for(auto& t : mThreads) t.join();
        }

        uint64_t addTask(LPFunc f, void* args) {
            Task t;
            t.func = f;
            t.args = args;
            std::lock_guard<std::mutex> lk(mMutex);
            uint64_t sz = mTasks.size();
            mTasks.push(std::move(t));
            mCondActive.notify_all();
            return sz;
        }
        uint64_t addTaskUnsafe(LPFunc f, void* args) {
            Task t;
            t.func = f;
            t.args = args;
            uint64_t sz = mTasks.size();
            mTasks.push(std::move(t));
            return sz;
        }

    private:
        uint64_t mNumOfThreads;
        std::vector<std::thread> mThreads;
        
        std::condition_variable mCondActive;
        std::condition_variable mCondMaster;
        std::mutex mMutex;
        uint64_t mWaiting;
        bool mPoolStopped, mPoolRunning;

        std::queue<Task> mTasks;

        std::atomic<uint32_t> mHelperId;

        bool inline canProceed() {
            return ((mPoolRunning && !mTasks.empty()) || mPoolStopped);
        }

        void worker(uint32_t tid) {
            //cerr << tid << endl;
            while(true) {
                std::unique_lock<std::mutex> lk(mMutex);
                ++mWaiting;
                //cerr << ">" << endl;
                // signal for synchronization!!! - all threads are waiting
                if (mWaiting == mNumOfThreads && mTasks.empty()) mCondMaster.notify_all();
                if (!canProceed()) mCondActive.wait(lk, [tid,this]{return canProceed();});
                --mWaiting;
                if (mPoolStopped) { lk.unlock(); return; }
                
                auto cTask = mTasks.front();
                mTasks.pop();
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
                std::unique_lock<std::mutex> lk(mMutex);
                //std::cerr << ">" << std::endl;
                // all threads are waiting and no jobs are pending
                if (mWaiting == mNumOfThreads ||  mTasks.empty()) {
                    lk.unlock();
                    return;
                }
                if (mPoolStopped || !mPoolRunning) { lk.unlock(); return; }

                auto cTask = mTasks.front();
                mTasks.pop();
                lk.unlock();

                //std::cerr << "got job!" << std::endl;

                // run the function passing the thread ID
                cTask.func(mNumOfThreads, tid, cTask.args);
            }
        }

};

#endif