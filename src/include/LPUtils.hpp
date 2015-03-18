#ifndef __LP_UTILS__
#define __LP_UTILS__

#pragma once

#include <chrono>
#include <thread>
#include <vector>
#include <algorithm>

#include "agner/vectorclass.h"
#include "aligned_allocator.hpp"

namespace lp {

#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)

#define ALWAYS_INLINE  __attribute__((always_inline)) inline

#define lp_EQUAL(x, y) (!((x)^(y)))
#define lp_EQUAL2(x, y) (!(1 + ~(x) + (y)))

namespace utils {

    template<typename T>
    bool ALWAYS_INLINE exists(T *a, const size_t sz, T val) {
        return std::find_if(a, a+sz, [&](const T& c) { return c == val; }) != a+sz;
    }

    template<typename T>
    uint64_t find_zero(T *a, const size_t sz) {
        for (size_t i = 0; i<sz; ++i) if (a[i] == 0) return i;
        return sz;
    }

}

namespace simd {   
    typedef uint64_t auint __attribute__ ((__aligned__(16)));
    template<unsigned int asz = 16> using auint64 __attribute__ ((__aligned__(asz)))= uint64_t;
    
    template<typename T, unsigned int asz = 16> using a16_t __attribute__ ((__aligned__(asz)))= T;
    
    uint64_t ALWAYS_INLINE or_all(auint64<16> *__restrict__ a, const unsigned int sz) {
        auint64<16> ored = 0;
        for (unsigned int i=0; i<sz; ++i) {
            ored |= a[i];
        }
        return ored;
    }
    uint64_t ALWAYS_INLINE or_all(std::vector<uint64_t> a) {
        auint ored = 0; const unsigned int sz = a.size();
        for (unsigned int i=0; i<sz; ++i) {
            ored |= a[i];
        }
        return ored;
    }
    
    template<typename A>
    bool ALWAYS_INLINE exists(std::vector<A> const& a, const A val) {
        auint mask = 0;
        const unsigned int sz = a.size();
        for (unsigned int i=0; i<sz; ++i) {
            mask |= a[i] == val;
        }
        return mask;
    }
    /* 
    template<typename T>
    bool ALWAYS_INLINE exists_early(a16_t<T> *__restrict__ a, const size_t sz, T val) {
        for (size_t i=0; i<sz; ++i) {
            if (a[i] == val) return true;
        }
        return false;
    }
    */
    template<typename T>
    bool ALWAYS_INLINE exists(a16_t<T> *__restrict__ a, const size_t sz, T val) {
        return std::find_if(a, a+sz, [&](const T& c) { return c == val; }) != a+sz;
    }
    
    bool ALWAYS_INLINE exists(a16_t<uint64_t> *__restrict__ a, const size_t sz, uint64_t val) {
        Vec2uq veca;
        if (sz&1) {
            for (size_t i=1; i<sz; i += 2) {
                veca.load(a+i);
                if (horizontal_or(veca == val)) return true;
            }
            return a[0] == val;
        } else {
            for (size_t i=0; i<sz; i += 2) {
                veca.load(a+i);
                if (horizontal_or(veca == val)) return true;
            }
            return false;
        }
    }
    bool ALWAYS_INLINE exists_avx(a16_t<uint64_t> *__restrict__ a, const size_t sz, uint64_t val) {
        Vec4uq veca;
        switch (sz&3) {
            case 0:
                for (size_t i=0; i<sz; i += 4) {
                    veca.load(a+i);
                    if (horizontal_or(veca == val)) return true;
                }
                return false;
            case 1:
                for (size_t i=1; i<sz; i += 4) {
                    veca.load(a+i);
                    if (horizontal_or(veca == val)) return true;
                }
                return a[0] == val;
            case 2:
                for (size_t i=2; i<sz; i += 4) {
                    veca.load(a+i);
                    if (horizontal_or(veca == val)) return true;
                }
                return a[0] == val || a[1] == val;
            case 3:
                for (size_t i=3; i<sz; i += 4) {
                    veca.load(a+i);
                    if (horizontal_or(veca == val)) return true;
                }
                return a[0] == val || a[1] == val || a[2] == val;
        }
        return false;
    }
    
    // SPECIFIC functions for the aligned allocator types I use
    template<typename T>
    using aligned_vector = std::vector<T, aligned_allocator<T, 16>>;
    
    template<typename A>
    bool ALWAYS_INLINE exists(aligned_vector<A> const& a, const A val) {
        return exists(a.data(), a.size(), val);
    }
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
