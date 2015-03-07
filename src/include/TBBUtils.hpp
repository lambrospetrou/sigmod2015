#ifndef __LP_TBB_UTILS__
#define __LP_TBB_UTILS__

#pragma once

#include <tbb/tbb.h>
#include <tbb/parallel_for.h>

//https://software.intel.com/en-us/blogs/2007/11/08/have-a-fish-how-break-from-a-parallel-loop-in-tbb

//! Like a blocked_range, but becomes immediately empty if "stop" flag is true.
template<typename Value>
class cancelable_range {
    tbb::blocked_range<Value> my_range;
    volatile bool& my_stop;

    public:
    // Constructor for client code
    /** Range becomes empty if stop==true. */
    cancelable_range( int begin, int end, int grainsize, volatile bool& stop ) :
        my_range(begin,end,grainsize),
        my_stop(stop) {}

    //! Splitting constructor used by parallel_for
    cancelable_range( cancelable_range& r, tbb::split ) :
        my_range(r.my_range,tbb::split()),
        my_stop(r.my_stop) {}

    //! Cancel the range.
    void cancel() const {my_stop=true;}

    //! True if range is empty.
    /** Range is empty if there is request to cancel the range. */
    bool empty() const {return my_stop || my_range.empty();}

    //! True if range is divisible
    /** Range becomes indivisible if there is request to cancel the range. */
    bool is_divisible() const {return !my_stop && my_range.is_divisible();}

    //! Initial value in range.
    Value begin() const {return my_range.begin();}

    //!One past last value in range
    // Note that end()==begin() if there is request to cancel the range.
    // The value of end() may change asynchronously if another thread cancels the range.
    Value end() const {return my_stop ? my_range.begin() : my_range.end();}

};

#endif
