#ifndef __LP_BOUNDEDSRSWQUEUE__
#define __LP_BOUNDEDSRSWQUEUE__

#include <vector>
#include <mutex>
#include <condition_variable>
#include <cstdint>

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
        mCondFull.wait(lk, [this]{return !isFull();});
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
        mCondFull.notify_one();    
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
/*
    inline bool isEmpty() const { 
        if (mEnqIndex < mDeqIndex) return mEnqIndex + mCapacity == mDeqIndex; 
        else return mEnqIndex == mDeqIndex; 
    }
    inline bool isFull() const {
        if (mEnqIndex < mDeqIndex) return (mEnqIndex + mCapacity) - mDeqIndex == mMaxSize;
        else return mEnqIndex - mDeqIndex == mMaxSize;
    }
*/
};

#endif
