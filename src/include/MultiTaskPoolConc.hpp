#ifndef __LP_MULTITASKPOOL_CONC__
#define __LP_MULTITASKPOOL_CONC__

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

#include "tbb/tbb.h"

class MultiTaskPoolConc {
    private:

        //template<typename T> using LPFunc = std::function<void(uint32_t, uint32_t, T)>;
        typedef std::function<void(uint32_t, uint32_t, void*)> LPFunc;

        struct Task {
            LPFunc func;
            void* args;
        };

    public:
        MultiTaskPoolConc(uint32_t tsz) : mNumOfThreads(tsz), 
            mWaiting(0), mPoolStopped(false), mPoolRunning(false), mHelperId(tsz) {
        }

        void initThreads() {
            for (uint32_t i=0; i<mNumOfThreads; ++i) {
                mThreads.push_back(std::thread(&MultiTaskPoolConc::worker, this, i));
            }
        }

        inline void startAll() {
            std::unique_lock<std::mutex> lk(mMutex);
            mPoolRunning = true;
            lk.unlock();
            mCondActive.notify_all();
        }

        inline void waitAllAndStop() {
            std::unique_lock<std::mutex> lk(mMutex);
            mCondMaster.wait(lk, [this]{
                    return (mWaiting == mNumOfThreads && mTasks.empty());
                    });
            mPoolRunning = false;
            lk.unlock();
        }
        inline void waitAll() {
            //std::unique_lock<std::mutex> lk(mMutex);
            //mCondMaster.wait(lk, [this]{
            //        return (mWaiting == mNumOfThreads && mTasks.empty());
            //        });
            //lk.unlock();
            while (!mTasks.empty() || mWaiting != mNumOfThreads) std::this_thread::yield();
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
            mTasks.push(std::move(t));
            //std::lock_guard<std::mutex> lk(mMutex);
            mCondActive.notify_one();
        }
        inline void addTaskUnsafe(LPFunc f, void* args) {
            Task t;
            t.func = f;
            t.args = args;
            mTasks.push(std::move(t));
        }

    private:
        uint64_t mNumOfThreads;
        std::vector<std::thread> mThreads;
        
        std::condition_variable mCondActive;
        std::condition_variable mCondMaster;
        std::mutex mMutex;
        uint64_t mWaiting;
        bool mPoolStopped, mPoolRunning;

        tbb::concurrent_queue<Task> mTasks;

        std::atomic<uint32_t> mHelperId;

        bool inline canProceed() {
            return ((mPoolRunning && !mTasks.empty()) || mPoolStopped);
        }

        void worker(uint32_t tid) {
            //cerr << tid << endl;
            while(true) {
                //std::unique_lock<std::mutex> lk(mMutex);
                ++mWaiting;
                //cerr << ">" << endl;
                // signal for synchronization!!! - all threads are waiting
                if (mWaiting == mNumOfThreads && mTasks.empty()) mCondMaster.notify_one();
                if (!canProceed()) { 
                    std::unique_lock<std::mutex> lk(mMutex);
                    mCondActive.wait(lk, [tid,this]{return canProceed();});
                    lk.unlock();
                }
                --mWaiting;
                if (mPoolStopped) { /*lk.unlock();*/ return; }
                
                Task cTask; bool popped = false;
                while (!(popped=mTasks.try_pop(cTask)) && !mTasks.empty()); 
                if (!popped || mTasks.empty()) continue;
                //lk.unlock();
                //cerr << "-" << endl;

                // run the function passing the thread ID
                //cTask.func(mNumOfThreads, tid, std::move(cTask.args));
                cTask.func(mNumOfThreads, tid, cTask.args);
                // wait after the execution for the other threads
                //cerr << "2" << endl;
            }
        }

};

#endif
