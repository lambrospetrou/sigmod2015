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
        }

        // returns the bucket that will hold the tuples for thie given transaction
        // to make the insertions faster
        ALWAYS_INLINE Bucket* bucketNext(uint64_t trid) {
            if (unlikely(mBucketSize - mBuckets.back().trsize == 0)) {
                mBuckets.emplace_back(trid, trid, 1);
                //mBuckets.push_back({});
                //Bucket& lb = mBuckets.back();
                //lb.trsize = 1; lb.trmin = trid; lb.trmax = trid;
                return &mBuckets.back();
            } else {
                //Bucket& lb = mBuckets.back();
                //++lb.trsize; lb.trmax = trid;
                //if (lb.trsize==1) lb.trmin = trid;
                return mBuckets.back().setMax(trid);
            }
        }

        void ALWAYS_INLINE sortAll() {
            for (Bucket& b : mBuckets) { b.sortByVal(); }
        }

        // sorts the buckets that contain all transactions from trfrom and greater
        void ALWAYS_INLINE sortFrom(uint64_t trfrom) {
            auto bend = mBuckets.end();
            auto bfrom = std::lower_bound(mBuckets.begin(), bend, trfrom, BTRLess);
            while (bfrom<bend) {
                (bfrom++)->sortByVal(); 
            }
        }

        // returns the buckets that contain all the transactions from trfrom to trto
        // result [first, last)
        std::pair<Bucket*, Bucket*> buckets(uint64_t trfrom, uint64_t trto) {
            auto be = mBuckets.data()+mBuckets.size();
            //auto trf = std::lower_bound(mBuckets.data(), be, trfrom, BTRLess);
            //return {trf, std::upper_bound(trf, be, trto, BTRLess)};
            auto trf = mBuckets.data();
            for (;(be-trf>0) & (trf->trmax < trfrom);) ++trf;
            auto trt = trf;
            for (;(be-trt>0) & (trt->trmin <= trto);) ++trt;
            return {trf, trt};
        }

    private:

        std::vector<Bucket> mBuckets; // there should be at least one bucket ALWAYS
        //uint64_t mLBT; // least bucket transaction = the Bucket.from of the first bucket

};

#endif
