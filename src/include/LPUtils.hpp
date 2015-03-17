#ifndef __LP_UTILS__
#define __LP_UTILS__

#pragma once

#include <chrono>
#include <thread>
#include <vector>

namespace lp {

#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)

#define ALWAYS_INLINE  __attribute__((always_inline)) inline

#define lp_EQUAL(x, y) (!((x)^(y)))
#define lp_EQUAL2(x, y) (!(1 + ~(x) + (y)))

    
    typedef uint64_t auint __attribute__ ((__aligned__(16)));
    template<unsigned int asz = 16> using auint64 __attribute__ ((__aligned__(asz)))= uint64_t;
    
    uint64_t ALWAYS_INLINE lp_OR(auint64<16> *__restrict__ a, const unsigned int sz) {
        auint64<16> ored = 0;
        for (unsigned int i=0; i<sz; ++i) {
            ored |= a[i];
        }
        return ored;
    }
    uint64_t ALWAYS_INLINE lp_OR(std::vector<uint64_t> a) {
        auint ored = 0; const unsigned int sz = a.size();
        for (unsigned int i=0; i<sz; ++i) {
            ored |= a[i];
        }
        return ored;
    }

    inline void lp_spin_sleep(std::function<bool ()> pred) {
        do { std::this_thread::yield(); } while (!pred());
    }
    inline void lp_spin_sleep(std::chrono::microseconds us = std::chrono::microseconds(0)) {
        if (us == std::chrono::microseconds(0)) {
            std::this_thread::yield();
        } else {
            std::this_thread::sleep_for(us);
        }
    }
    inline void lp_spin_for(std::chrono::microseconds us = std::chrono::microseconds(0)) {
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
