#ifndef __LP_LPTIMER__
#define __LP_LPTIMER__

#pragma once

#include <iostream>
#include <sys/time.h>

struct LPTimer_t {
    uint64_t validations;
    uint64_t validationsProcessing;
    uint64_t validationsProcessingIndex;
    uint64_t transactions;
    uint64_t transactionsIndex;
    uint64_t updateIndex;
    uint64_t queryIndex;
    
    uint64_t flushes;
    uint64_t forgets;
    uint64_t reading;
    uint64_t readingTotal;
    LPTimer_t() : validations(0), validationsProcessing(0), validationsProcessingIndex(0), transactions(0), transactionsIndex(0), updateIndex(0), 
                   queryIndex(0), flushes(0), forgets(0), reading(0), readingTotal(0) {}
    uint64_t getChronoMicro() {      
        // return nanoseconds
        struct timeval start;
        gettimeofday(&start, NULL);
        // tv_sec = seconds | tv_usecs = microseconds
        return (start.tv_sec * 1000000LL) + start.tv_usec;
    }
    uint64_t getChronoMicro(uint64_t start) {      
        return getChronoMicro() - start;
    }
    uint64_t getChrono() {
        return getChronoMicro();
    }
    uint64_t getChrono(uint64_t start) {
        return getChronoMicro() - start;
    } 
}; 

std::ostream& operator<< (std::ostream& os, const LPTimer_t& t) {
    os << "LPTimer [val: " << t.validations << " val-proc: " << t.validationsProcessing << " == proc: " << t.validationsProcessingIndex << " trans: " << t.transactions << " trans-index: " << t.transactionsIndex << " upd-index: " << t.updateIndex << " q-index: " << t.queryIndex << " flush: " << t.flushes << " forget: " << t.forgets << " reads: " << t.reading << "/" << t.readingTotal<< "]" << std::endl;
    return os;
}


#endif
