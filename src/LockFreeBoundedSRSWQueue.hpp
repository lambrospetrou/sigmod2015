#ifndef __LP_LF_BOUNDEDSRSWQUEUE__
#define __LP_LF_BOUNDEDSRSWQUEUE__

#include <vector>
#include <mutex>
#include <condition_variable>
#include <cstdint>
#include <chrono>
#include <thread>
#include <functional>
#include <atomic>

// previous:: Only for Single Reader - Single Writer
template<class T>
class LockFreeBoundedSRSWQueue {
public:
    LockFreeBoundedSRSWQueue(uint64_t sz) : mCurSize(0), mEnqIndex(0), mDeqIndex(0), mMaxSize(sz), mCapacity(sz+1),
        mQ(std::vector<T>(mCapacity)) {}

    T& reqNextEnq() { 
        // Wait until the buffer is at least half empty
        lp_spin_sleep([this]{return !isHalfFull();});
        return mQ[mEnqIndex]; 
    }
    void registerEnq() {
        // might not be necessary
        mEnqIndex = (mEnqIndex+1) % mCapacity;
        ++mCurSize; // atomically
    }
    
    T& reqNextDeq() {
        lp_spin_sleep([this]{return !isEmpty();});
        return mQ[mDeqIndex];
    }
    void registerDeq() {
        mDeqIndex = (mDeqIndex + 1) % mCapacity;
        --mCurSize;
    }

private:
    std::atomic<uint64_t> mCurSize;  // this is just for speedup in isFull and isEmpty
    uint64_t mEnqIndex; // points to where the next empty slot is
    uint64_t mDeqIndex; // points to the next available message
    uint64_t mMaxSize;  // the maximum number iof messages in the queue
    uint64_t mCapacity; // the total storage
    std::vector<T> mQ;        

    inline bool isEmpty() const { return mCurSize.load() == 0; }
    inline bool isFull() const { return mCurSize.load() == mMaxSize; }
    inline bool isHalfFull() const { return mCurSize.load() >= (mMaxSize>>1); }
    
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
