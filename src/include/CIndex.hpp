#ifndef __LP_COLUMN_INDEX__
#define __LP_COLUMN_INDEX__

#include "LPUtils.hpp"
#include "DoubleIter.hpp"
#include "aligned_allocator.hpp"

#include <cstdint>
#include <vector>
#include <algorithm>

#include <iostream>

class CIndex {

    template<typename T>
    using vector_a = std::vector<T, aligned_allocator<T, 16>>;
    //using vector_a = std::vector<T>;
    
    using tuple_t = uint64_t*;
    static constexpr size_t mBucketSize = 256;
    
    public:
        struct Meta_t {
            //uint32_t tpl_id;
            uint64_t trans_id;
            tuple_t tuple;
        };
        
        struct Bucket {
            vector_a<uint64_t> values; // might be btree maybe if faster than binary_search
            vector_a<Meta_t> meta;
            uint64_t trmin;
            uint64_t trmax;
            size_t trsize;
            // other possible statistics for these values goes here

            Bucket() : trmin(0), trmax(0), trsize(0) {
                values.reserve(mBucketSize); values.reserve(mBucketSize);
            }
            Bucket(uint64_t _min, uint64_t _max, uint64_t sz) : trmin(_min), trmax(_max), trsize(sz) {
                values.reserve(mBucketSize); values.reserve(mBucketSize);
            }

            Meta_t* rawMeta() { return meta.data(); }
            uint64_t* rawValues() { return values.data(); }

            void ALWAYS_INLINE notifyInsertBatch(size_t sz) { values.reserve(values.size()+sz); meta.reserve(meta.size()+sz); }

            void ALWAYS_INLINE insert(uint64_t trid, tuple_t tpl, uint64_t val) {
                values.push_back(val);
                meta.push_back({trid, tpl});
            }

            ALWAYS_INLINE Bucket* setMax(uint64_t trid) {
                trmin = (trsize > 0) ? trmin : trid;
                ++trsize; trmax = trid;
                return this;
            }

            void sortByVal() {
                std::sort(SIter<uint64_t, Meta_t>(values.data(), meta.data()), 
                    SIter<uint64_t, Meta_t>(values.data()+values.size(), meta.data()+values.size()));
            }

            std::pair<Meta_t*, Meta_t*> equal_range(uint64_t v) {
                //auto vb = values.data(), ve = values.data()+values.size();
                //auto rp = std::equal_range(vb, ve, v);
                //auto mb = meta.data();
                //auto vb = values.begin(), ve = values.end();
                //auto rp = std::equal_range(vb, ve, v);
                auto vb = values.data();
                auto rp = std::equal_range(vb, vb+values.size(), v);
                
                //for (auto vv : values) std::cerr << vv << " ";
                //std::cerr << "\nval: " << v << " sz: " << values.size() << " res: " << (rp.second-rp.first) << std::endl;
                return {meta.data()+(rp.first-vb), meta.data()+(rp.second-vb)};
            }
        };
  
        struct BTRLess_t {
            ALWAYS_INLINE bool operator() (const Bucket& left, const Bucket& right) {
                return left.trmax < right.trmax;
            }
            ALWAYS_INLINE bool operator() (const Bucket& o, uint64_t target) {
                return o.trmax < target;
            }
            ALWAYS_INLINE bool operator() (uint64_t target, const Bucket& o) {
                return target < o.trmin;
            }
        } BTRLess;

        CIndex() {
            mBuckets.resize(1);
            mBB = mBuckets.data();
            mBE = mBB+1;
        }

        // returns the bucket that will hold the tuples for thie given transaction
        // to make the insertions faster
        ALWAYS_INLINE Bucket* bucketNext(uint64_t trid) {
            if (unlikely(mBucketSize - mBuckets.back().trsize == 0)) {
                mBuckets.emplace_back(trid, trid, 1);
                mBB = mBuckets.data(); mBE = mBuckets.data()+mBuckets.size();
                return &mBuckets.back();
            } else {
                return mBuckets.back().setMax(trid);
            }
        }

        void ALWAYS_INLINE sortAll() {
            for (Bucket& b : mBuckets) { b.sortByVal(); }
        }

       
        /////////////////////
        /////////
        ///////// UTILITY FUNCTIONS 
        /////////
        /////////////////////
        ALWAYS_INLINE Bucket* lp_lower_bound(uint64_t trid) {
            auto trt = std::lower_bound(mBB, mBE, trid, BTRLess);
            //auto trt = mBB;
            //for (;(mBE-trt>0) & (trt->trmax < trid); ++trt);
            return trt;
        }
        template<typename Iter>
        ALWAYS_INLINE Bucket* lp_upper_bound(uint64_t trid, Iter bb) {
            auto trt = std::upper_bound(bb, mBE, trid, BTRLess);
            //auto trt = bb;
            //for (;(mBE-trt>0) & (trt->trmin <= trid); ++trt);
            return trt;
        }
        /////////////////////

        // sorts the buckets that contain all transactions from trfrom and greater
        void ALWAYS_INLINE sortFrom(uint64_t trfrom) {
            auto trt = lp_lower_bound(trfrom);
            while (trt<mBE) {
                (trt++)->sortByVal(); 
            }
        }

        // returns the buckets that contain all the transactions from trfrom to trto
        // result [first, last)
        std::pair<Bucket*, Bucket*> buckets(uint64_t trfrom, uint64_t trto) {
            //////
            //auto be = mBuckets.data()+mBuckets.size();
            //auto trf = mBB;
            //for (;(mBE-trf>0) & (trf->trmax < trfrom); ++trf);
            //auto trt = trf;
            //for (;(mBE-trt>0) & (trt->trmin <= trto); ++trt);
            //return {trf, trt};
            auto trf = lp_lower_bound(trfrom);
            return {trf, lp_upper_bound(trto, trf)};
        }

        ALWAYS_INLINE std::pair<Bucket*, Bucket*> buckets() {
            return {mBB, mBE};
        }

        ALWAYS_INLINE void erase(uint64_t trto) {
            auto trt = lp_lower_bound(trto);
            if (trt->trmax - trto == 0) {
                //std::cerr << "d";
                mBuckets.erase(mBuckets.begin(), mBuckets.begin()+(trt-mBB));
                mBB = mBuckets.data(); mBE = mBuckets.data()+mBuckets.size();
            } else if (trt > mBB) {
                //std::cerr << "d";
                mBuckets.erase(mBuckets.begin(), mBuckets.begin()+(trt-1-mBB));    
                mBB = mBuckets.data(); mBE = mBuckets.data()+mBuckets.size();
            }
        }

    private:

        std::vector<Bucket> mBuckets; // there should be at least one bucket ALWAYS
        Bucket *mBB; // mBuckets.data()
        Bucket *mBE; // mBuckets.data()+mBuckets.size()
};

#endif
