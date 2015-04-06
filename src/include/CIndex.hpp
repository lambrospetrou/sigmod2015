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
    //using vector_a = std::vector<T, aligned_allocator<T, 16>>;
    using vector_a = std::vector<T>;
    
    using tuple_t = uint64_t*;
    static constexpr size_t BUCKET_TUPLES_LIMIT = ((size_t)1)<<11;
    static constexpr size_t BUCKET_TRANS_LIMIT = 128;
    static constexpr size_t BUCKET_PRIMARY_LIMIT = 512;
    
    public:
        struct Meta_t {
            uint64_t value;
            uint64_t trans_id;
            tuple_t tuple;
        };
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
        };
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
        };
        struct MTrLess_t {
            ALWAYS_INLINE bool operator() (const Meta_t& left, const Meta_t& right) {
                return left.trans_id < right.trans_id;
            }
            ALWAYS_INLINE bool operator() (const Meta_t& o, uint64_t target) {
                return o.trans_id < target;
            }
            ALWAYS_INLINE bool operator() (uint64_t target, const Meta_t& o) {
                return target < o.trans_id;
            }
        };
        
        struct Bucket {
            //vector_a<uint64_t> values; // might be btree maybe if faster than binary_search
            vector_a<Meta_t> meta;
            uint64_t trmin;
            uint64_t trmax;
            size_t trsize;
            // other possible statistics for these values goes here

            Bucket() : trmin(0), trmax(0), trsize(0) {
                //values.reserve(mBucketSize); 
                meta.reserve(BUCKET_TUPLES_LIMIT);
            }
            Bucket(uint64_t _min, uint64_t _max, uint64_t sz) : trmin(_min), trmax(_max), trsize(sz) {
                //values.reserve(mBucketSize); meta.reserve(mBucketSize);
                meta.reserve(BUCKET_TUPLES_LIMIT);
            }

            //void ALWAYS_INLINE notifyInsertBatch(size_t sz) { values.reserve(values.size()+sz); meta.reserve(meta.size()+sz); }
           
            // return a pointer to the next empty position of those allocated now
            /*
            std::pair<uint64_t*, Meta_t*> ALWAYS_INLINE resizeAndGetPtr(size_t sz) { 
                size_t oldsz = values.size();
                values.resize(values.size()+sz); 
                meta.resize(meta.size()+sz); 
                return {values.data()+oldsz, meta.data()+oldsz};
            }
            */
            std::pair<uint64_t*, Meta_t*> ALWAYS_INLINE resizeAndGetPtr(size_t sz) { 
                size_t oldsz = meta.size();
                meta.resize(oldsz+sz); 
                return {nullptr, meta.data()+oldsz};
            }

            void ALWAYS_INLINE insert(uint64_t trid, tuple_t tpl, uint64_t val) {
                //values.push_back(val);
                //meta.push_back({trid, tpl});
                meta.push_back({val, trid, tpl});
            }

            ALWAYS_INLINE Bucket* setMax(uint64_t trid) {
                trmin = (trsize > 0) ? trmin : trid;
                ++trsize; trmax = trid;
                return this;
            }

            void sortByVal() {
                //std::sort(SIter<uint64_t, Meta_t>(values.data(), meta.data()), 
                //    SIter<uint64_t, Meta_t>(values.data()+values.size(), meta.data()+values.size()));
                std::sort(meta.data(), meta.data()+meta.size(), MLess_t());
            }
            void sortByValTrans() {
                //std::sort(SIter<uint64_t, Meta_t>(values.data(), meta.data()), 
                //    SIter<uint64_t, Meta_t>(values.data()+values.size(), meta.data()+values.size()));
                std::sort(meta.data(), meta.data()+meta.size(), MVTLess_t());
            }

            std::pair<Meta_t*, Meta_t*> tuples() {
                return {meta.data(), meta.data()+meta.size()};
            }

            std::pair<Meta_t*, Meta_t*> equal_range(uint64_t v, uint64_t trfrom, uint64_t trto) {
                auto trp = std::equal_range(meta.data(), meta.data()+meta.size(), v, MVTLess_t());
                auto ub = std::upper_bound(trp.first, trp.second, trto, MTrLess_t());
                auto lb = std::lower_bound(trp.first, ub, trfrom, MTrLess_t());
                return {lb, ub};
            }
            std::pair<Meta_t*, Meta_t*> equal_range(uint64_t v, uint64_t trid) {
                auto trp = std::equal_range(meta.data(), meta.data()+meta.size(), v, MVTLess_t());
                return {trp.first, std::upper_bound(trp.first, trp.second, trid, MTrLess_t())};
            }
            std::pair<Meta_t*, Meta_t*> equal_range(uint64_t v) {
                //auto vb = values.data();
                //auto rp = std::equal_range(vb, vb+values.size(), v);
                //return {meta.data()+(rp.first-vb), meta.data()+(rp.second-vb)};
                return std::equal_range(meta.data(), meta.data()+meta.size(), v, MVTLess_t());
            }
            std::pair<Meta_t*, Meta_t*> lower_bound(uint64_t v) {
                //auto vb = values.data();
                //auto lb = std::lower_bound(vb, vb+values.size(), v);
                //return {meta.data()+(lb-vb), meta.data()+meta.size()};
                auto mb = meta.data(), me = mb + meta.size(); 
                return {std::lower_bound(mb, me, v, MVTLess_t()), me};
            }
            std::pair<Meta_t*, Meta_t*> lower_bound_left(uint64_t v) {
                //auto vb = values.data();
                //auto lb = std::lower_bound(vb, vb+values.size(), v);
                //return {meta.data(), meta.data()+(lb-vb)};
                auto mb = meta.data(), me = mb + meta.size(); 
                return {mb, std::lower_bound(mb, me, v, MVTLess_t())};
            }
            std::pair<Meta_t*, Meta_t*> upper_bound(uint64_t v) {
                //auto vb = values.data();
                //auto ub = std::upper_bound(vb, vb+values.size(), v);
                //return {meta.data()+(ub-vb), meta.data()+meta.size()};
                auto mb = meta.data(), me = mb + meta.size(); 
                return {std::upper_bound(mb, me, v, MVTLess_t()), me};
            }
            std::pair<Meta_t*, Meta_t*> upper_bound_left(uint64_t v) {
                //auto vb = values.data();
                //auto ub = std::upper_bound(vb, vb+values.size(), v);
                //return {meta.data(), meta.data()+(ub-vb)};
                auto mb = meta.data(), me = mb + meta.size(); 
                return {mb, std::upper_bound(mb, me, v, MVTLess_t())};
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

        // takes into account only the number of transactions
        ALWAYS_INLINE bool is_newbucket_primary() {
            return (unlikely(mBuckets.empty() || (mBuckets.back().trsize >= BUCKET_PRIMARY_LIMIT)));
        }
        // takes into account the number of tuples too
        ALWAYS_INLINE bool is_newbucket() { 
            return (mBuckets.empty() || 
                    (mBuckets.back().trsize >= BUCKET_TRANS_LIMIT) ||
                    (mBuckets.back().meta.size() >= BUCKET_TUPLES_LIMIT));
        }
        // returns the bucket that will hold the tuples for thie given transaction
        // to make the insertions faster
        ALWAYS_INLINE Bucket* bucketNext(uint64_t trid, bool isPrimary = false) {
            //if (unlikely(mBuckets.empty() || (mBucketSize - mBuckets.back().trsize == 0))) {
            if (unlikely(isPrimary ? is_newbucket_primary() : is_newbucket())) {
                mBuckets.emplace_back(trid, trid, 1);
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

#define BINARY_SZ 16

        ALWAYS_INLINE Bucket* lp_lower_bound(uint64_t trid) {
            if (unlikely(mBuckets.size() > BINARY_SZ)) {
                return std::lower_bound(BB(), BE(), trid, BTRLess);
            } else {
                auto mBB = BB(), mBE = BE();
                auto trt = mBB;
                for (;(mBE-trt>0) && (trt->trmax < trid); ++trt);
                return trt;
            }
        }
        template<typename Iter>
        ALWAYS_INLINE Bucket* lp_upper_bound(uint64_t trid, Iter bb) {
            if (unlikely(mBuckets.size() > BINARY_SZ)) {
                return std::upper_bound(bb, BE(), trid, BTRLess);
            } else {
                auto mBE = BE();
                auto trt = bb;
                for (;(mBE-trt>0) && (trt->trmin <= trid); ++trt);
                return trt;
            }
        }
        ALWAYS_INLINE Bucket* lp_upper_bound(uint64_t trid) {
            return lp_upper_bound(trid, BB());
        }
        /////////////////////

        // sorts the buckets that contain all transactions from trfrom and greater
        void ALWAYS_INLINE sortFrom(uint64_t trfrom, bool noTransSort = false) {
            auto mBE = BE();
            auto trt = lp_lower_bound(trfrom);
            if (likely(!noTransSort)) {
                while (trt<mBE) {
                    (trt++)->sortByValTrans(); 
                }
            } else {
                while (trt<mBE) {
                    (trt++)->sortByVal(); 
                }
            }
        }

        // returns the buckets that contain all the transactions from trfrom to trto
        // result [first, last)
        std::pair<Bucket*, Bucket*> buckets(uint64_t trfrom, uint64_t trto) {
            auto trf = lp_lower_bound(trfrom);
            return {trf, lp_upper_bound(trto, trf)};
        }

        ALWAYS_INLINE std::pair<Bucket*, Bucket*> buckets() {
            return {BB(), BE()};
        }

        ALWAYS_INLINE void erase(uint64_t trto) {
            if (mBuckets.size() < 3) return;
            //if (unlikely(mBuckets.size() == 1)) return;
            auto trt = lp_lower_bound(trto);
            if (trt->trmax - trto == 0) {
                //std::cerr << "d";
                mBuckets.erase(mBuckets.begin(), mBuckets.begin()+(trt-BB()));
            } else if (trt > BB()) {
                //std::cerr << "d";
                mBuckets.erase(mBuckets.begin(), mBuckets.begin()+(trt-1-BB()));    
            }
        }

        size_t size() const { return mBuckets.size(); }

    private:
        std::vector<Bucket> mBuckets; 
};

#endif
