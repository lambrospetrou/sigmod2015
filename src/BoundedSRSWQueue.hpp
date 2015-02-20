#ifndef __LP_BOUNDEDSRSWQUEUE__
#define __LP_BOUNDEDSRSWQUEUE__

#include <vector>
#include <mutex>
#include <condition_variable>
#include <cstdint>
#include <chrono>
#include <thread>
#include <functional>

// NOW IT IS MRMW !!!
// previous:: Only for Single Reader - Single Writer
template<class T>
class BoundedSRSWQueue {
public:
    BoundedSRSWQueue(uint64_t sz) : mCurSize(0), mEnqIndex(0), mDeqIndex(0), mMaxSize(sz), mCapacity(sz+1),
        mQ(std::vector<T>(mCapacity)) {}

    T& reqNextEnq() { 
        // Wait until main() sends data
        std::unique_lock<std::mutex> lk(mMutex);
        mCondFull.wait(lk, [this]{return !isHalfFull();});
        lk.unlock();
        return mQ[mEnqIndex]; 
    }
    void registerEnq() {
        // might not be necessary
        std::lock_guard<std::mutex> lk(mMutex);
        mEnqIndex = (mEnqIndex+1) % mCapacity;
        ++mCurSize;
        mCondEmpty.notify_one();
    }
    
    T& reqNextDeq() {
        std::unique_lock<std::mutex> lk(mMutex);
        mCondEmpty.wait(lk, [this]{return !isEmpty();});
        lk.unlock();
        return mQ[mDeqIndex];
    }
    void registerDeq() {
        std::lock_guard<std::mutex> lk(mMutex);
        mDeqIndex = (mDeqIndex + 1) % mCapacity;
        --mCurSize;
        if (!isHalfFull()) mCondFull.notify_one();    
    }

private:
    uint64_t mCurSize;  // this is just for speedup in isFull and isEmpty
    uint64_t mEnqIndex; // points to where the next empty slot is
    uint64_t mDeqIndex; // points to the next available message
    uint64_t mMaxSize;  // the maximum number iof messages in the queue
    uint64_t mCapacity; // the total storage
    std::vector<T> mQ;        

    //std::mutex mMutexEnq;
    //std::mutex mMutexDeq;
    std::mutex mMutex;
    std::condition_variable mCondFull;
    std::condition_variable mCondEmpty;

    inline bool isEmpty() const { return mCurSize == 0; }
    inline bool isFull() const { return mCurSize == mMaxSize; }
    inline bool isHalfFull() const { return mCurSize >= (mMaxSize>>1); }
    
    // busy waiting with yield in order to allow other threads waiting in execution
    // queue to be executed
    inline void lp_spin_sleep(std::function<bool ()> pred) {
        do { std::this_thread::yield(); } while (!pred());
    }
    inline void lp_spin_sleep(std::chrono::microseconds us = std::chrono::microseconds(0)) {
        if (us == std::chrono::microseconds(0)) {
            std::this_thread::yield();
        } else {
            auto start = std::chrono::high_resolution_clock::now();
            auto tend = start + us;
            do { std::this_thread::yield(); }
            while (std::chrono::high_resolution_clock::now() < tend);
        }
    }

};

#endif
