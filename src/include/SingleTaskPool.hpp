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
        SingleTaskPool(uint32_t tsz, void(*func)(uint32_t, uint32_t)) : mNumOfThreads(tsz), mThreadsActivity(0), 
            mWaiting(0), mPoolStopped(false), mPoolRunning(false) {
            mMaskAll = (static_cast<uint64_t>(1) << mNumOfThreads) - 1;
            mFunc = func;
        }

        void initThreads() {
            for (uint32_t i=0; i<mNumOfThreads; ++i) {
                mThreads.push_back(std::thread(&SingleTaskPool::worker, this, i));
            }
        }

        void startAll() {
            std::unique_lock<std::mutex> lk(mMutex);
            mPoolRunning = true;
            lk.unlock();
            mCondActive.notify_all();
        }
        void startAll(void (*func)(uint32_t, uint32_t)) {
            std::unique_lock<std::mutex> lk(mMutex);
            mFunc = func;
            mPoolRunning = true;
            lk.unlock();
            mCondActive.notify_all();
        }

        void waitAll() {
            std::unique_lock<std::mutex> lk(mMutex);
            mCondMaster.wait(lk, [this]{
                    return (mThreadsActivity == mMaskAll) && (mWaiting == mNumOfThreads);
                    });
            mPoolRunning = false;
            mThreadsActivity = 0;
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
        uint64_t mMaskAll;
        uint64_t mNumOfThreads;
        uint64_t mThreadsActivity;
        std::vector<std::thread> mThreads;
        
        std::condition_variable mCondActive;
        std::condition_variable mCondMaster;
        std::mutex mMutex;
        uint64_t mWaiting;
        bool mPoolStopped, mPoolRunning;

        // the function to be called
        void (*mFunc)(uint32_t, uint32_t);

        void worker(uint32_t tid) {
            uint64_t tMask = static_cast<uint64_t>(1)<<tid;
            //cerr << tid << endl;
            while(true) {
                std::unique_lock<std::mutex> lk(mMutex);
                ++mWaiting;
                //cerr << ">" << endl;
                // signal for synchronization!!! - all threads are waiting
                if (mWaiting == mNumOfThreads) mCondMaster.notify_all();
                mCondActive.wait(lk, [tid,tMask,this]{return ((mThreadsActivity & tMask)==0 && mPoolRunning) || mPoolStopped;});
                --mWaiting;
                mThreadsActivity |= tMask;
                lk.unlock();
                // check if the pool is still active
                if (mPoolStopped) { 
                    //cerr << "e" <<tid << endl; 
                    return; 
                }
                //cerr << "-" << endl;

                // run the function passing the thread ID
                mFunc(mNumOfThreads, tid);
                // wait after the execution for the other threads
                //cerr << "2" << endl;
            }
        }

};

#endif
