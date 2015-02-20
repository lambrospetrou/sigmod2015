#ifndef __LP_SINGLETASKPOOL__
#define __LP_SINGLETASKPOOL__

#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <cstdlib>
#include <cstdint>
#include <cstring>

class SingleTaskPool {
    public:
        SingleTaskPool(uint32_t tsz) : mNumOfThreads(tsz), mPoolStopped(false), mPoolRunning(false), mWaiting(0) {
        }

        void initThreads(void (*func)(uint32_t, uint32_t)) {
            mThreadsActivity.resize(mNumOfThreads);
            memset(mThreadsActivity.data(), 0, sizeof(uint32_t)*mNumOfThreads);
            for (uint32_t i=0; i<mNumOfThreads; ++i) {
                mThreads.push_back(std::thread(&SingleTaskPool::worker, this, i, func));
            }
        }

        void startAll() {
            std::unique_lock<std::mutex> lk(mMutex);
            mPoolRunning = true;
            lk.unlock();
            mCondActive.notify_all();
        }

        void waitAll() {
            std::unique_lock<std::mutex> lk(mMutex);
            mCondMaster.wait(lk, [this]{
                    for (auto ta : mThreadsActivity) if (ta == 0) return false;
                    return mWaiting == mNumOfThreads;
                    });
            mPoolRunning = false;
            memset(mThreadsActivity.data(), 0, sizeof(uint32_t)*mNumOfThreads);
            lk.unlock();
        }

        void destroy() {
            std::unique_lock<std::mutex> lk(mMutex);
            mPoolStopped = true;
            lk.unlock();
            mCondActive.notify_all();
            mCondMaster.notify_all();
            for(auto& t : mThreads) t.join();
        }

    private:
        uint32_t mNumOfThreads;
        std::vector<std::thread> mThreads;
        std::vector<uint32_t> mThreadsActivity;
        bool mPoolStopped, mPoolRunning;
        std::condition_variable mCondActive;
        std::mutex mMutex;
        uint64_t mWaiting;
        std::condition_variable mCondMaster;

        void worker(uint32_t tid, void (*func)(uint32_t, uint32_t)) {
            //cerr << tid << endl;
            while(true) {
                std::unique_lock<std::mutex> lk(mMutex);
                ++mWaiting;
                //cerr << ">" << endl;
                // signal for synchronization!!! - all threads are waiting
                if (mWaiting == mNumOfThreads) mCondMaster.notify_all();
                mCondActive.wait(lk, [tid, this]{return ((mThreadsActivity[tid]==0 && mPoolRunning) || mPoolStopped);});
                --mWaiting;
                lk.unlock();
                mThreadsActivity[tid] = 1;
                // check if the pool is still active
                if (mPoolStopped) { 
                    //cerr << "e" <<tid << endl; 
                    return; 
                }
                //cerr << "-" << endl;

                // run the function passing the thread ID
                func(mNumOfThreads, tid);
                // wait after the execution for the other threads
                //cerr << "2" << endl;
            }
        }
};

#endif
