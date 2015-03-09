#pragma once

#include <algorithm>
#include <iostream>

template<typename IterA, typename IterB>
class SIter {
    public:

        SIter(IterA aIter, IterB bIter) : aIter(aIter), bIter(bIter) {}

        // we move to next position is just moving both:
        SIter& operator ++() {
            ++aIter; ++bIter;
            return *this;
        }
        SIter operator ++(int) {
            SIter rv = *this;
            ++aIter; ++bIter;
            return rv;
        }
        SIter& operator --() {
            --aIter; --bIter;
            return *this;
        }
        SIter operator --(int) {
            SIter rv = *this;
            --aIter; --bIter;
            return rv;
        }
        SIter operator + (std::ptrdiff_t cc) const
        {
            SIter rv = *this;
            rv.aIter += cc;
            rv.bIter += cc;
            return rv;
        }
        SIter operator - (std::ptrdiff_t cc) const
        {
            SIter rv = *this;
            rv.aIter -= cc;
            rv.bIter -= cc;
            return rv;
        }
        std::ptrdiff_t operator - (SIter other) const
        {
            return aIter - other.aIter;
        }
        struct value_type {
            typename IterA::value_type a; typename IterB::value_type b;
            bool operator < (const value_type& other) const {
                return a < other.a; // this is the place where way of sorting is defined!!!!
            }
        };
        struct reference {
            IterA a;
            IterB b;
            reference& operator = (const value_type& o) 
            {
                *a = o.a;
                *b = o.b;
                return *this;
            }
            operator value_type() const {
                value_type rv = { *a, *b };
                return rv;
            }
            reference& operator = (const reference& other)
            {
                *a = *other.a;
                *b = *other.b;
                return *this;
            }
            bool operator < (const reference& other) const {
                return *a < *other.a; 
            }
            bool operator < (const value_type& other) const {
                return *a < other.a; 
            }
        };

        reference operator * () {
            reference rv = { aIter, bIter };
            return rv;
        }
        bool operator == (const SIter& other) const
        {
            return aIter == other.aIter; // don't need to compare bIter - shall be in sync
        }
        bool operator != (const SIter& other) const
        {
            return aIter != other.aIter; // don't need to compare bIter - shall be in sync
        }
        bool operator < (const SIter& other) const
        {
            return aIter < other.aIter; // don't need to compare bIter - shall be in sync
        }
        // I bet you don't need pointer operator -> ()
        typedef std::random_access_iterator_tag iterator_category;
        typedef std::ptrdiff_t difference_type;
        typedef reference pointer;


    private:
        IterA aIter; 
        IterB bIter;
};
