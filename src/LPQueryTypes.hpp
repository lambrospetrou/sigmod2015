#ifndef __LP_TYPES__
#define __LP_TYPES__

#include "ReferenceTypes.hpp"
#include <vector>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <ostream>
#include <iostream>
#include <utility>
#include <algorithm>

namespace lp {

    using namespace std;

    struct LPQuery {
        /// The relation
        uint32_t relationId;
        /// The number of bound columns
        uint32_t columnCount;
        /// The bindings
        std::vector<Query::Column> predicates;
        /// States whether this query can ever be true - will be set to false later by our filtering
        bool satisfiable;

        LPQuery() : relationId(-1), columnCount(0), predicates(std::vector<Query::Column>()), satisfiable(true) {}
        LPQuery(const Query& q) : relationId(q.relationId), columnCount(q.columnCount), predicates(std::vector<Query::Column>(q.columnCount)), satisfiable(true) {
            memcpy(predicates.data(), q.columns, sizeof(Query::Column)*columnCount);
            std::sort(predicates.begin(), predicates.end());
            predicates.resize(std::distance(predicates.begin(), std::unique(predicates.begin(), predicates.end())));
            //if (columns.size() != columnCount) cerr << "diff: " << columnCount-columns.size() << endl;
            columnCount = predicates.size();
        }
        // this is the default - operator< of Column
        static bool QCSortCol (const Query::Column& left, const Query::Column& right) {
            if (left.column < right.column) return true;
            else if (right.column < left.column) return false;
            else if (left.op < right.op) return true;
            else if (right.op < left.op) return false;
            else return left.value < right.value;    
        }
        static bool QCSortOp (const Query::Column& left, const Query::Column& right) {
            if (left.op < right.op) return true;
            else if (right.op < left.op) return false;
            else if (left.column < right.column) return true;
            else if (right.column < left.column) return false;
            else return left.value < right.value;    
        }
        static bool LPQuerySizeLessThan(const LPQuery& left, const LPQuery& right) {
            //if (left.satisfiable && !right.satisfiable) return true;
            //else if (right.satisfiable < !left.satisfiable) return false;
            //else return (left.columnCount < right.columnCount);
            return (left.columnCount < right.columnCount);
        }
    };

    bool operator< (const LPQuery& left, const LPQuery& right) {
        if (left.relationId < right.relationId) return true;
        else if (right.relationId < left.relationId) return false;
        else if (left.columnCount < right.columnCount) return true;
        else return left.predicates < right.predicates;
    }

    bool operator== (const LPQuery& left, const LPQuery& right)  {
        if (left.relationId != right.relationId) return false;
        if (left.columnCount != right.columnCount) return false;
        return left.predicates == right.predicates;
    }
    

    
    bool operator< (const Query::Column& left, const Query::Column& right) {
        if (left.column < right.column) return true;
        else if (right.column < left.column) return false;
        else if (left.op < right.op) return true;
        else if (right.op < left.op) return false;
        else return left.value < right.value;    
    }
    
    bool operator== (const Query::Column& left, const Query::Column& right) {
        if (left.column != right.column) return false;
        else if (left.op != right.op) return false;
        else return left.value == right.value;    
    }
    
    struct LPValidation {
        uint64_t validationId;
        uint64_t from,to;
        std::vector<LPQuery> queries;
        LPValidation(const ValidationQueries& v, std::vector<LPQuery> q)
            : validationId(v.validationId), from(v.from), to(v.to), queries(move(q)) {}
        LPValidation(uint64_t vid, uint64_t fr, uint64_t t, std::vector<LPQuery> q)
            : validationId(vid), from(fr), to(t), queries(q) {}
    };
    
    // the following will be used for the query filtering and quick rejection
    namespace validation {

        typedef Query::Column::Op Op;
        typedef Query::Column Column;

        struct Satisfiability {
            bool pastOps[6];
            uint64_t eq=UINT64_MAX, lt = UINT64_MAX, leq = UINT64_MAX, gt = 0, geq = 0;
            std::vector<uint64_t> neq; 
            Satisfiability():eq(UINT64_MAX),lt(UINT64_MAX), leq(UINT64_MAX), gt(0), geq(0) {
                memset(pastOps, 0, 6);
            }
            void reset() {
                eq=UINT64_MAX; lt = UINT64_MAX; leq = UINT64_MAX; gt = 0; geq = 0;
                memset(pastOps, 0, 6);
                neq.clear();
            }
        };

        bool isQueryColumnUnsolvable(Column& p, Satisfiability& sat) {
            switch (p.op) {
                case Op::Equal:
                    // already found an equality check
                    if (sat.pastOps[Op::Equal]) return true;
                    sat.pastOps[Op::Equal] = true; 
                    sat.eq = p.value;
                    break;
                case Op::NotEqual:
                    // both equality and inequality with same value
                    if (sat.pastOps[Op::Equal] && sat.eq == p.value) return true;
                    sat.pastOps[Op::NotEqual] = true;
                    //sat.neq.push_back(p.value);
                    break;
                case Op::Less:
                    if (sat.pastOps[Op::Equal] && sat.eq >= p.value) return true;
                    sat.pastOps[Op::Less] = true;
                    if (p.value < sat.lt) { sat.lt = p.value; sat.leq = p.value - 1; }
                    break;
                case Op::LessOrEqual:
                    if (sat.pastOps[Op::Equal] && sat.eq > p.value) return true;
                    sat.pastOps[Op::LessOrEqual] = true;
                    if (p.value < sat.leq) { sat.leq = p.value; sat.lt = p.value + 1; }
                    break;
                case Op::Greater:
                    if (sat.pastOps[Op::Equal] && sat.eq <= p.value) return true;
                    sat.pastOps[Op::Greater] = true;
                    if (p.value > sat.gt) { sat.gt = p.value; sat.geq = p.value + 1; }
                    break;
                case Op::GreaterOrEqual:
                    if (sat.pastOps[Op::Equal] && sat.eq < p.value) return true;
                    sat.pastOps[Op::GreaterOrEqual] = true;
                    if (p.value > sat.geq) { sat.geq = p.value; sat.gt = p.value - 1; }
                    break;
            }

            // check for equality and constrasting ranges
            if (sat.pastOps[Op::Equal] && (sat.eq < sat.gt || sat.eq > sat.lt))
                return true;

            // check non-overlapping ranges
            if (sat.lt - sat.gt <= 0) return true;

            return false;
        }

        bool isQueryUnsolvable(LPQuery& q) {
            if (q.predicates.empty()) return false;
            Satisfiability sat;
            uint64_t lastCol = UINT64_MAX;
            for (auto& p : q.predicates) {
                if (p.column != lastCol) { sat.reset(); lastCol=p.column; }
                if (isQueryColumnUnsolvable(p, sat)) return true;
            }
            return false;
        }

        bool unsolvable(LPValidation& v) {
            // NOTE:: I ASSUME THAT THE PREDICATES ARE UNIQUE IN EACH LPQuery
            // NOTE:: I ASSUME THAT PREDICATES ARE SORTED ON THEIR COLUMNS FIRST
            uint32_t unsatisfied = 0;
            for (auto& q : v.queries) {
                if (isQueryUnsolvable(q)) { 
                    q.satisfiable = false; ++unsatisfied; 
                }
            }
            if (v.queries.size() == unsatisfied) return true;
            // since we have some queries that are candidates for confliction
            // we will sort them to be in front of the queries and be able to break
            /*uint64_t left=0, right=v.queries.size()-1, lastSwapPos=v.queries.size()-1;
              for (; left < right; ) {
              for(; v.queries[left].satisfiable && left<right; ) ++left;
              for(; !v.queries[right].satisfiable && left<right; ) --right;
              if (left < right) {
              std::swap(v.queries[left], v.queries[right]);
              lastSwapPos = right;
              ++left; --right;
              }
              }
              (void)lastSwapPos;*/
            return false;
        }
    }
}

#endif
