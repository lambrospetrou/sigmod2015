#ifndef __LP_TYPES__
#define __LP_TYPES__

#include "ReferenceTypes.hpp"
#include "LPUtils.hpp"
#include <vector>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <ostream>
#include <iostream>
#include <utility>
#include <algorithm>
#include <memory>

#include "asm/asmlib.h"

namespace lp {

    typedef Query::Column::Op Op;
    typedef Query::Column Column;

    struct ColumnCompColOnly_t {
        inline bool operator() (const Query::Column& left, const Query::Column& right) {
            return left.column < right.column;    
        }
    } ColumnCompColOnly;
    struct ColumnCompCol_t {
        inline bool operator() (const Query::Column& left, const Query::Column& right) {
            if (left.column < right.column) return true;
            else if (right.column < left.column) return false;
            else if (left.op < right.op) return true;
            else if (right.op < left.op) return false;
            else return left.value < right.value;    
        }
    } ColumnCompCol;
    struct ColumnCompColEq_t {
        inline bool operator() (const Query::Column& left, const Query::Column& right) {
            if (left.column != right.column) return false;
            else if (left.op != right.op) return false;
            else return left.value == right.value;    
        }
    } ColumnCompColEq;
    struct ColumnCompOp_t {
        inline bool operator()(const Query::Column *left, const Query::Column *right) {
            if (left->op < right->op) return true;
            else if (right->op < left->op) return false;
            else if (left->column < right->column) return true;
            else if (right->column < left->column) return false;
            else return left->value < right->value;    
        }
        inline bool operator()(const Query::Column& left, const Query::Column& right) {
            if (left.op < right.op) return true;
            else if (right.op < left.op) return false;
            else if (left.column < right.column) return true;
            else if (right.column < left.column) return false;
            else return left.value < right.value;    
        }
    } ColumnCompOp;

    struct ColumnCompQuality_t {
        inline bool operator() (const Query::Column& left, const Query::Column& right) {
            if (left.op == Op::Equal && right.op != Op::Equal) return true;
            else if (left.op != Op::Equal && right.op == Op::Equal) return false;
            else if (left.op == Op::NotEqual && right.op != Op::NotEqual) return false;
            else if (left.op != Op::NotEqual && right.op == Op::NotEqual) return true;
            return (left.column < right.column);
        }
    } ColumnCompQuality;

    struct LPQuery {
        /// The relation
        uint32_t relationId;
        /// The number of bound columns
        uint32_t columnCount;
        /// The raw query as read by the input
        Query *rawQuery;
        uint32_t colCountUniq;
        //std::vector<Query::Column> predicates;
        /// States whether this query can ever be true - will be set to false later by our filtering
        //bool satisfiable;

        LPQuery() : relationId(-1), columnCount(0), rawQuery(nullptr), colCountUniq(0) {}
        explicit LPQuery(Query *q) : relationId(q->relationId), columnCount(q->columnCount), rawQuery(q), colCountUniq(q->columnCount) {
            //if (q->columnCount == 0) return;
            //std::sort(q->columns, q->columns+q->columnCount, ColumnCompCol);
            //auto colEnd = std::unique(q->columns, q->columns+q->columnCount, ColumnCompColEq);
            //colCountUniq = std::distance(q->columns, colEnd);
        }
    };

    struct LPQueryCompSize_t {
        inline bool operator()(const LPQuery& left, const LPQuery& right) {
            return (left.columnCount < right.columnCount);
        }
    } LPQueryCompSizeLess;
    struct LPQueryCompSize2_t {
        inline bool operator()(const LPQuery& left, const LPQuery& right) {
            return (left.colCountUniq < right.colCountUniq);
        }
    } LPQueryCompUniqSize;
/*
    struct LPQueryCompQuality_t {
        inline bool operator()(const LPQuery& left, const LPQuery& right) {
            if (left.columnCount == 0 && right.columnCount > 0) return true;
            else if (left.columnCount > 0 && right.columnCount == 0) return false;
            return ColumnCompQuality(left.predicates[0], right.predicates[0]);
        }
    } LPQueryCompQuality;
*/
    struct ReceivedMessage; // forward declaration - declared in the ReaderIO header

    struct LPValidation {
        std::vector<LPQuery> queries;
        ReceivedMessage *rawMsg;
        
        uint64_t validationId;
        uint64_t from,to;
        
        LPValidation() {}
        LPValidation(uint64_t vid, uint64_t fr, uint64_t t, ReceivedMessage* msg, std::vector<LPQuery> q)
            : queries(move(q)), rawMsg(msg), validationId(vid), from(fr), to(t) {}

        ~LPValidation() {}
    };


    namespace query {

        // return the new number of valid predicates
        inline uint32_t __attribute__((always_inline)) preprocess(Query& rq) {
            if (unlikely(rq.columnCount == 0)) return 0;
            std::sort(rq.columns, rq.columns+rq.columnCount, ColumnCompCol);
            return std::distance(rq.columns, std::unique(rq.columns, rq.columns+rq.columnCount, ColumnCompColEq));
        }
        inline uint32_t __attribute__((always_inline)) preprocess(LPQuery& lpq) {
            auto q = lpq.rawQuery;
            std::sort(q->columns, q->columns+q->columnCount, ColumnCompCol);
            auto colEnd = std::unique(q->columns, q->columns+q->columnCount, ColumnCompColEq);
            lpq.colCountUniq = std::distance(q->columns, colEnd);
            return lpq.colCountUniq;
        }

        struct Satisfiability {
            bool pastOps[6];
            uint64_t eq=UINT64_MAX, lt = UINT64_MAX, leq = UINT64_MAX, gt = 0, geq = 0;
            Satisfiability():eq(UINT64_MAX),lt(UINT64_MAX), leq(UINT64_MAX), gt(0), geq(0) {
                A_memset(pastOps, 0, 6);
            }
            inline void reset() {
                eq=UINT64_MAX; lt = UINT64_MAX; leq = UINT64_MAX; gt = 0; geq = 0;
                A_memset(pastOps, 0, 6);
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
                    //sat.pastOps[Op::NotEqual] = true; // TODO - check what to do
                    break;
                case Op::Less:
                    if (sat.pastOps[Op::Equal] && sat.eq >= p.value) return true;
                    sat.pastOps[Op::Less] = true;
                    if (p.value < sat.lt) { sat.lt = p.value; sat.leq = p.value - 1; }
                    //else { p.value = sat.lt; }
                    break;
                case Op::LessOrEqual:
                    if (sat.pastOps[Op::Equal] && sat.eq > p.value) return true;
                    sat.pastOps[Op::LessOrEqual] = true;
                    if (p.value < sat.leq) { sat.leq = p.value; sat.lt = p.value + 1; }
                    //else { p.value = sat.leq; }
                    break;
                case Op::Greater:
                    if (sat.pastOps[Op::Equal] && sat.eq <= p.value) return true;
                    sat.pastOps[Op::Greater] = true;
                    if (p.value > sat.gt) { sat.gt = p.value; sat.geq = p.value + 1; }
                    //else { p.value = sat.gt; }
                    break;
                case Op::GreaterOrEqual:
                    if (sat.pastOps[Op::Equal] && sat.eq < p.value) return true;
                    sat.pastOps[Op::GreaterOrEqual] = true;
                    if (p.value > sat.geq) { sat.geq = p.value; sat.gt = p.value - 1; }
                    //else { p.value = sat.geq; }
                    break;
            }

            // check for equality and constrasting ranges
            if (sat.pastOps[Op::Equal] && (sat.eq < sat.gt || sat.eq > sat.lt))
                return true;

            // check non-overlapping ranges
            if (sat.lt - sat.gt <= 0) return true;

            return false;
        }

        bool isQueryUnsolvable(Column *colBegin, Column *colEnd) {
            if (colBegin == colEnd) return false;
            Satisfiability sat;
            uint64_t lastCol = UINT64_MAX;
            for (; colBegin != colEnd; ++colBegin) {
                if (colBegin->column != lastCol) { sat.reset(); lastCol=colBegin->column; }
                if (isQueryColumnUnsolvable(*colBegin, sat)) return true;
            }
            return false;
        }

        bool inline satisfiable(Query* q, uint32_t colCountUniq) {
            if (0 == colCountUniq) return true;
            Column *qc = const_cast<Column*>(q->columns);
            auto colBegin = qc, colEnd = qc + colCountUniq; //colEnd = qc + q->columnCount;
            if (isQueryUnsolvable(colBegin, colEnd)) {
                // the query is not-satisfiable so it should be skipped-pruned
                return false;
            }
            //std::partial_sort(colBegin, colBegin+std::min(q->columnCount, (uint32_t)2), colEnd, ColumnCompQuality);
            std::partial_sort(colBegin, colBegin+std::min(colCountUniq, (uint32_t)2), colEnd, ColumnCompQuality);
            return true;
        }
        bool inline satisfiable(LPQuery& q) { 
            /*
            Column * qc = const_cast<Column*>(q.rawQuery->columns);
            auto colBegin = qc, colEnd = qc + q.colCountUniq;
            if (isQueryUnsolvable(colBegin, colEnd)) {
                // the query is not-satisfiable so it should be skipped-pruned
                return false;
            }
            std::partial_sort(colBegin, colBegin+std::min(q.colCountUniq, (uint32_t)4), colEnd, ColumnCompQuality);
            return true;
            */
            return satisfiable(q.rawQuery, q.colCountUniq);
        }
/*
        bool parse(const Query *q, uint32_t relCols, LPQuery *nQ, uint32_t tid = 0) {
            (void)relCols; (void)tid;
            Column * qc = const_cast<Column*>(q->columns);
            
            //   std::cerr << std::endl;
            //   for (auto c=q->columns, cLimit=c+q->columnCount; c!=cLimit; ++c) {
            //   std::cerr << c->column << ":" << c->op << ":" << c->value << " ";
               }
             
            //std::vector<Column> preds(q->columns, q->columns+q->columnCount);

            // sort the columns by column first in order to remove uniques and check satisfiability
            std::sort(qc, qc+q->columnCount, ColumnCompCol);
            auto colEnd = std::unique(qc, qc+q->columnCount, ColumnCompColEq);
            auto colBegin = qc;
            uint64_t uniqSz = std::distance(colBegin, colEnd);

            if (isQueryUnsolvable(colBegin, colEnd)) {
                // the query is not-satisfiable so it should be skipped-pruned
                nQ->columnCount = 0;
                return false;
            }

            // the query is satisfiable so insert it into the predicates of the passed in LPQuery
            std::partial_sort(colBegin, colBegin+std::min(uniqSz, (uint64_t)2), colEnd, ColumnCompQuality);
            //nQ->predicates.insert(nQ->predicates.begin(), colBegin, colEnd);
            nQ->predicates.reserve(uniqSz);
            nQ->predicates.resize(uniqSz);
            memcpy(nQ->predicates.data(),colBegin, sizeof(Column)*uniqSz);
            nQ->columnCount = uniqSz;
            //std::cerr << "unique: " << uend-q->columns << std::endl; 

            //for (auto c: preds) std::cerr << c.column << ":" << c.op << ":" << c.value << " ";
            return true;
        }
*/
        


    }

    // the following will be used for the query filtering and quick rejection
    namespace validation {
        const static uint32_t MSK32_L16 = 0x0000ffff;
        const static uint32_t MSK32_H16 = 0xffff0000;

        uint32_t inline __attribute__((always_inline)) packRelCol(uint32_t rel, uint32_t col) {
            return ((rel & MSK32_L16) << 16) | (col & MSK32_L16);
        }
        void inline __attribute__((always_inline)) unpackRelCol(uint32_t pck, uint32_t& rel, uint32_t& col) {
            rel = ((pck & MSK32_H16) >> 16);
            col = (pck & MSK32_L16);
        }
        uint32_t inline __attribute__((always_inline)) unpackRel(uint32_t pck) {
            return ((pck & MSK32_H16) >> 16);
        }
    }
}

#endif
