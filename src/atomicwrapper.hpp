#ifndef __LP_ATOMIC_WRAPPER__
#define __LP_ATOMIC_WRAPPER__

#include <atomic>

template <typename T>
class atomic_wrapper {
    private:
        std::atomic<T> _a;
        //char _padding[64-sizeof(std::atomic<T>)];
    public:
        atomic_wrapper()
            :_a() {}

        atomic_wrapper(const std::atomic<T> &a)
            :_a(a.load()){}

        atomic_wrapper(const atomic_wrapper &other)
            :_a(other._a.load()){}

        atomic_wrapper &operator=(const atomic_wrapper &other) {
            _a.store(other._a.load());
            return *this;
        }

        void store(T desired) {_a.store(desired);}
        void store(T desired) volatile {_a.store(desired);}
        T load() const { return _a.load(); }
        T load() const volatile { return _a.load(); }

        std::atomic<T>& operator=(T desired) { _a.store(desired); return _a; }
        std::atomic<T>& operator=(T desired) volatile { _a.store(desired); return _a; }

        operator bool() const { return _a.load(); }
};

#endif
