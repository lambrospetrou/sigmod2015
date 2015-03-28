#pragma once

#include <algorithm>
#include <iostream>

template<typename TA, typename TB>
class SIter {
    public:
        SIter(TA* aIter, TB* bIter) : aIter(aIter), bIter(bIter) {}

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
        SIter& operator += (std::ptrdiff_t cc) {
            aIter += cc;
            bIter += cc;
            return *this;
        }
        SIter operator - (std::ptrdiff_t cc) const
        {
            SIter rv = *this;
            rv.aIter -= cc;
            rv.bIter -= cc;
            return rv;
        }
        SIter& operator -= (std::ptrdiff_t cc) {
            aIter -= cc;
            bIter -= cc;
            return *this;
        }
        std::ptrdiff_t operator - (SIter other) const
        {
            return aIter - other.aIter;
        }
        struct value_type {
            TA a; TB b;
            bool operator < (const value_type& other) const {
                return a < other.a; // this is the place where way of sorting is defined!!!!
            }
        };
        struct reference {
            friend void swap(reference a, reference b) {
                using std::swap;
                swap(*a.a, *b.a);
                swap(*a.b, *b.b);
            }
            TA* a;
            TB* b;
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
        TA* aIter; 
        TB* bIter;
};




/*
int main() {
    int a[10] = {10,9,8,7,6,5,4,3,2,1};
    double b[10] = {1,2,3,4,5,6,7,8,9,10};

    SIter beginIter(a, b);
    SIter endIter(a + 10, b + 10);

    std::sort(beginIter, endIter);
    for (int i = 0; i < 10; ++i) {
        std::cout << a[i] << "->" << b[i] << "\n";
    }
}
*/

