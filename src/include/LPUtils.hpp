#ifndef __LP_UTILS__
#define __LP_UTILS__

#pragma once

#include <chrono>
#include <thread>

namespace lp {

#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)

#define ALWAYS_INLINE  __attribute__((always_inline)) inline

#define lp_EQUAL(x, y) (!((x)^(y)))
#define lp_EQUAL2(x, y) (!(1 + ~(x) + (y)))


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
