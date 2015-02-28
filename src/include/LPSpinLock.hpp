#ifndef __LP_SPINLOCK__
#define __LP_SPINLOCK__

#include <atomic>

class LPSpinLock {
    private:
        std::atomic_flag mState = ATOMIC_FLAG_INIT;
    public:
        LPSpinLock() {}
        void inline lock() {
            while (mState.test_and_set(std::memory_order_acquire));
        }
        void inline unlock() {
            mState.clear(std::memory_order_release);
        }
};
/*
class LPSpinLock {
    private:
        enum SpinLockState { Locked, Unlocked };
        std::atomic<SpinLockState> mState;
    public:
        LPSpinLock() : mState(Unlocked) {}
        void inline lock() {
            while (mState.exchange(Locked, std::memory_order_acquire) == Locked);
        }
        void inline unlock() {
            mState.store(Unlocked, std::memory_order_release);
        }
};
*/
class LPSpinLockGuard {
    LPSpinLock &mLock;
    public:
    explicit LPSpinLockGuard(LPSpinLock &lock) : mLock(lock) { mLock.lock(); }
    ~LPSpinLockGuard() { mLock.unlock(); }
};


#endif
