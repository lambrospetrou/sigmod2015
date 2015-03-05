#ifndef __LP_UTILS__
#define __LP_UTILS__

#pragma once

#include <chrono>
#include <thread>

namespace lp {

#define _likely(x)       __builtin_expect(!!(x), 1)
#define _unlikely(x)     __builtin_expect(!!(x), 0)

    inline bool __attribute__((always_inline)) unlikely(bool cond) {
        return (__builtin_expect(cond, 0));
    }
    inline bool __attribute__((always_inline)) likely(bool cond) {
        return (__builtin_expect(cond, 1));
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
