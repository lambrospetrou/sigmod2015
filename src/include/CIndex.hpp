#ifndef __LP_COLUMN_INDEX__
#define __LP_COLUMN_INDEX__

#include "LPUtils.hpp"
#include "DoubleIter.hpp"

#include <cstdint>
#include <vector>
#include <algorithm>

#include <iostream>

class CIndex {

    using tuple_t = uint64_t*;
    static constexpr size_t mBucketSize = 128;
    
    public:
        struct Meta_t {
            //uint32_t tpl_id;
            uint64_t trans_id;
            tuple_t tuple;
        };
        
        struct Bucket {
            std::vector<uint64_t> values; // might be btree maybe if faster than binary_search
            std::vector<Meta_t> meta;
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

            void insert(uint64_t trid, tuple_t tpl, uint64_t val) {
                values.push_back(val);
                meta.push_back({trid, tpl});
            }

            void sortByVal() {
                std::sort(SIter<uint64_t, Meta_t>(values.data(), meta.data()), 
                    SIter<uint64_t, Meta_t>(values.data()+values.size(), meta.data()+values.size()));
            }

            std::pair<Meta_t*, Meta_t*> equal_range(uint64_t v) {
                auto vb = values.data(), ve = values.data()+values.size();
                auto mb = meta.data();
                auto rp = std::equal_range(vb, ve, v);
                
                //for (auto vv : values) std::cerr << vv << " ";
                //std::cerr << "\nval: " << v << " sz: " << values.size() << " res: " << (rp.second-rp.first) << std::endl;
                return {mb+(rp.first-vb), mb+(rp.second-vb)};
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
        };

        CIndex() {
            mBuckets.resize(1);
        }

        // returns the bucket that will hold the tuples for thie given transaction
        // to make the insertions faster
        ALWAYS_INLINE Bucket* bucketNext(uint64_t trid) {
            if (mBuckets.back().trsize == mBucketSize) {
                mBuckets.emplace_back(trid, trid, 1);
                //mBuckets.push_back({});
                //Bucket& lb = mBuckets.back();
                //lb.trsize = 1; lb.trmin = trid; lb.trmax = trid;
            } else {
                Bucket& lb = mBuckets.back();
                if (lb.trsize == 0) lb.trmin = trid;
                ++lb.trsize; lb.trmax = trid;
            }
            return &mBuckets.back();
        }

        void ALWAYS_INLINE sortAll() {
            for (Bucket& b : mBuckets) { b.sortByVal(); }
        }

        // sorts the buckets that contain all transactions from trfrom and greater
        void ALWAYS_INLINE sortFrom(uint64_t trfrom) {
            auto bend = mBuckets.end();
            auto bfrom = std::lower_bound(mBuckets.begin(), bend, trfrom, BTRLess_t());
            while (bfrom<bend) {
                bfrom->sortByVal(); 
                ++bfrom;
            }
        }

        // returns the buckets that contain all the transactions from trfrom to trto
        std::pair<Bucket*, Bucket*> buckets(uint64_t trfrom, uint64_t trto) {
            auto bb = mBuckets.data();
            auto be = mBuckets.data()+mBuckets.size();
            auto trf = std::lower_bound(bb, be, trfrom, BTRLess_t());
            return {trf, std::upper_bound(trf, be, trto, BTRLess_t())};
        }

    private:

        std::vector<Bucket> mBuckets; // there should be at least one bucket ALWAYS
        //uint64_t mLBT; // least bucket transaction = the Bucket.from of the first bucket

};

#endif
