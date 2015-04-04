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
            //uint64_t value;
            uint64_t trans_id;
            tuple_t tuple;
        };
        /*
        struct MLess_t {
            ALWAYS_INLINE bool operator() (const Meta_t& left, const Meta_t& right) {
                return left.value < right.value;
            }
            ALWAYS_INLINE bool operator() (const Meta_t& o, uint64_t target) {
                return o.value < target;
            }
            ALWAYS_INLINE bool operator() (uint64_t target, const Meta_t& o) {
                return target < o.value;
            }
        } MetaLess;
        struct MVTLess_t {
            ALWAYS_INLINE bool operator() (const Meta_t& left, const Meta_t& right) {
                if (left.value < right.value) return true;
                else if (right.value < left.value) return false;
                else return left.trans_id < right.trans_id;
            }
            ALWAYS_INLINE bool operator() (const Meta_t& left, uint64_t target) {
                return left.value < target;
            }
            ALWAYS_INLINE bool operator() (uint64_t target, const Meta_t& right) {
                return target < right.value;
            }
        } MetaVTLess;
        */
        struct Bucket {
            vector_a<uint64_t> values; // might be btree maybe if faster than binary_search
            vector_a<Meta_t> meta;
            uint64_t trmin;
            uint64_t trmax;
            size_t trsize;
            // other possible statistics for these values goes here

            Bucket() : trmin(0), trmax(0), trsize(0) {
                values.reserve(mBucketSize); meta.reserve(mBucketSize);
            }
            Bucket(uint64_t _min, uint64_t _max, uint64_t sz) : trmin(_min), trmax(_max), trsize(sz) {
                values.reserve(mBucketSize); meta.reserve(mBucketSize);
            }

            void ALWAYS_INLINE notifyInsertBatch(size_t sz) { values.reserve(values.size()+sz); meta.reserve(meta.size()+sz); }

            void ALWAYS_INLINE insert(uint64_t trid, tuple_t tpl, uint64_t val) {
                values.push_back(val);
                meta.push_back({trid, tpl});
                //meta.push_back({val, trid, tpl});
            }

            ALWAYS_INLINE Bucket* setMax(uint64_t trid) {
                trmin = (trsize > 0) ? trmin : trid;
                ++trsize; trmax = trid;
                return this;
            }

            void sortByVal() {
                std::sort(SIter<uint64_t, Meta_t>(values.data(), meta.data()), 
                    SIter<uint64_t, Meta_t>(values.data()+values.size(), meta.data()+values.size()));
                //std::sort(meta.data(), meta.data()+meta.size(), MVTLess_t());
            }

            std::pair<Meta_t*, Meta_t*> equal_range(uint64_t v) {
                auto vb = values.data();
                auto rp = std::equal_range(vb, vb+values.size(), v);
                //for (auto vv : values) std::cerr << vv << " ";
                //std::cerr << "\nval: " << v << " sz: " << values.size() << " res: " << (rp.second-rp.first) << std::endl;
                return {meta.data()+(rp.first-vb), meta.data()+(rp.second-vb)};
                //return std::equal_range(meta.data(), meta.data()+meta.size(), v, MVTLess_t());
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
            if (unlikely(mBuckets.empty() || (mBucketSize - mBuckets.back().trsize == 0))) {
                mBuckets.emplace_back(trid, trid, 1);
                //mBB = mBuckets.data(); mBE = mBuckets.data()+mBuckets.size();
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

#define BB() (mBuckets.data())
#define BE() (mBuckets.data() + mBuckets.size())

        ALWAYS_INLINE Bucket* lp_lower_bound(uint64_t trid) {
            //return std::lower_bound(BB(), BE(), trid, BTRLess);
            auto mBB = BB(), mBE = BE();
            auto trt = mBB;
            for (;(mBE-trt>0) & (trt->trmax < trid); ++trt);
            return trt;
        }
        template<typename Iter>
        ALWAYS_INLINE Bucket* lp_upper_bound(uint64_t trid, Iter bb) {
            //return std::upper_bound(bb, BE(), trid, BTRLess);
            auto mBE = BE();
            auto trt = bb;
            for (;(mBE-trt>0) & (trt->trmin <= trid); ++trt);
            return trt;
        }
        ALWAYS_INLINE Bucket* lp_upper_bound(uint64_t trid) {
            return lp_upper_bound(trid, BB());
        }
        /////////////////////

        // sorts the buckets that contain all transactions from trfrom and greater
        void ALWAYS_INLINE sortFrom(uint64_t trfrom) {
            auto mBE = BE();
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
            return {BB(), BE()};
        }

        ALWAYS_INLINE void erase(uint64_t trto) {
            if (likely(mBuckets.size() < 4)) return;
            auto trt = lp_lower_bound(trto);
            if (trt->trmax - trto == 0) {
                //std::cerr << "d";
                mBuckets.erase(mBuckets.begin(), mBuckets.begin()+(trt-BB()));
                //mBuckets.erase(mBuckets.begin(), mBuckets.begin()+(trt-mBuckets.data()));
            } else if (trt > BB()) {
                //std::cerr << "d";
                mBuckets.erase(mBuckets.begin(), mBuckets.begin()+(trt-1-BB()));    
                //mBuckets.erase(mBuckets.begin(), mBuckets.begin()+(trt-1-mBuckets.data()));    
            }
        }

        size_t size() const { return mBuckets.size(); }

    private:

        std::vector<Bucket> mBuckets; // there should be at least one bucket ALWAYS
};

#endif
