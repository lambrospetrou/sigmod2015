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
            bool ALWAYS_INLINE exists_loop(T *a, const size_t sz, T val) {
                for (size_t i=0; i<sz; ++i) {
                    if (a[i] == val) return true;
                }
                return false;
            }
        template<typename T>
            bool ALWAYS_INLINE exists(T *a, const size_t sz, T val) {
                return std::find_if(a, a+sz, [&](const T& c) { return c == val; }) != a+sz;
            }

        template<typename T>
            uint64_t find_zero(T *a, const size_t sz) {
                for (size_t i = 0; i<sz; ++i) if (a[i] == 0) return i;
                return sz;
            }

        // https://schani.wordpress.com/2010/04/30/linear-vs-binary-search/
        bool ALWAYS_INLINE binary_cmov (const uint64_t *arr, const size_t sz, const uint64_t key) {
            size_t min = 0, max = sz;
            while (min < max) {
            //do {
                size_t middle = min + ((max-min) >> 1);
                asm ("cmpq %3, %2\n\tcmova %4, %0\n\tcmovbe %5, %1"
                        : "+r" (min),
                          "+r" (max)
                        : "r" (key), "g" (arr [middle]),
                        "g" (middle + 1), "g" (middle));
            
            }// while (min < max);
            return min < sz && lp_EQUAL(arr[min], key);
        }

    }

    namespace simd {   
        typedef uint64_t auint __attribute__ ((__aligned__(16)));
        template<unsigned int asz = 16> using auint64 __attribute__ ((__aligned__(asz)))= uint64_t;

        template<typename T, unsigned int asz = 16> using a16_t __attribute__ ((__aligned__(asz)))= T;

        void ALWAYS_INLINE and_left_auto(a16_t<uint8_t> *__restrict__ l, a16_t<uint8_t> *__restrict__ r, const size_t sz) {
            for (size_t i=0; i<sz; ++i) {
                l[i] &= r[i];
            }
        }
        void ALWAYS_INLINE and_left(a16_t<uint8_t> *__restrict__ l, a16_t<uint8_t> *__restrict__ r, const size_t sz) {
            const size_t extra = sz&3; // sz%4
            const size_t nsz = sz-extra;
            for (size_t i=0; i<nsz; i+=4) {
                l[i] &= r[i]; l[i+1] &= r[i+1]; l[i+2] &= r[i+2]; l[i+3] &= r[i+3];
            }
            switch (extra) {
                case 0: return;
                case 1: 
                    l[nsz] &= r[nsz]; return;
                case 2:
                    l += nsz; r += nsz;
                    l[0] &= r[0]; l[1] &= r[1]; return;
                case 3:
                    l += nsz; r += nsz;
                    l[0] &= r[0]; l[1] &= r[1]; l[2] &= r[2]; return;
            }
        }
        void ALWAYS_INLINE and_left_8(a16_t<uint8_t> *__restrict__ l, a16_t<uint8_t> *__restrict__ r, const size_t sz) {
            for (size_t i=0; i<sz; i+=8) {
                l[i] &= r[i];     l[i+1] &= r[i+1]; l[i+2] &= r[i+2]; l[i+3] &= r[i+3]; 
                l[i+4] &= r[i+4]; l[i+5] &= r[i+5]; l[i+6] &= r[i+6]; l[i+7] &= r[i+7];         
            }
        }
        void ALWAYS_INLINE and_left_4(a16_t<uint8_t> *__restrict__ l, a16_t<uint8_t> *__restrict__ r, const size_t sz) {
            for (size_t i=0; i<sz; i+=4) {
                l[i] &= r[i]; l[i+1] &= r[i+1]; l[i+2] &= r[i+2]; l[i+3] &= r[i+3]; 
            }
        }
        void ALWAYS_INLINE and_left_opt(a16_t<uint8_t> *__restrict__ l, a16_t<uint8_t> *__restrict__ r, const size_t sz) {
            Vec16uc veca, vecb;
            for (size_t i=0; i<sz; i+=16) {
                veca.load(l+i); vecb.load(r+i);
                veca &= vecb;
                veca.store(l+i);
            }
        }



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

        bool ALWAYS_INLINE exists(register const a16_t<uint64_t> *__restrict__ a, const size_t sz, const uint64_t val) {
            const size_t extra = sz&1; // sz%2
            //register a16_t<uint64_t> *a = (a16_t<uint64_t>*)__builtin_assume_aligned (_a, 16);
            register const a16_t<uint64_t> *aend = a+sz-extra;
            Vec2uq veca;
            for (; a<aend; ) {
                veca.load(a++); a++;
                if (horizontal_or(veca == val)) return true;
            }
            switch (extra) {
                case 0: 
                    return false;
                case 1: 
                    return a[0] == val;
            }
            return false;
        }
        bool ALWAYS_INLINE exists_avx(register const a16_t<uint64_t> *__restrict__ a, const size_t sz, const uint64_t val) {
            const size_t extra = sz&3; // sz%4
            //register a16_t<uint64_t> *a = (a16_t<uint64_t>*)__builtin_assume_aligned (_a, 16);
            register const a16_t<uint64_t> *aend = a+sz-extra;
            Vec4uq veca;
            for (; a<aend; a+=4) {
                veca.load(a); 
                if (horizontal_or(veca == val)) return true;
            }
            switch (extra) {
                case 0: 
                    return false;
                case 1: 
                    return *a == val;
                case 2:
                    return *a++ == val || *a == val;
                    //return *a == val || *(a+1) == val;
                case 3:
                    return *a++ == val || *a++ == val || *a == val;
                    //return horizontal_or(Vec4uq(*a, *(a+1), *(a+2), 0) == val);
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
