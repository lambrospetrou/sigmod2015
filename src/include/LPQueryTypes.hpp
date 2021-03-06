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
            if ((left.column ^ right.column)) return false;
            else if (left.op ^ right.op) return false;
            else return !(left.value ^ right.value);    
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
        ALWAYS_INLINE bool operator() (const Query::Column& left, const Query::Column& right) {
            /*
            if ((left.op == Op::Equal) & (right.op != Op::Equal)) return true;
            else if ((left.op != Op::Equal) & (right.op == Op::Equal)) return false;
            else if ((left.op == Op::NotEqual) & (right.op != Op::NotEqual)) return false;
            else if ((left.op != Op::NotEqual) & (right.op == Op::NotEqual)) return true;
            return (left.column < right.column);
            */
            /*
            if ((left.op == Op::Equal)) {
                if (right.op != Op::Equal) { return true; }
                else { return left.column < right.column; }
            } else {
                if (right.op == Op::Equal) { return false; } 
                else { return left.op > right.op; }
            }
            */
            if ((left.op) == (right.op)) { return left.column < right.column; }
            else if ((left.op == 0)) { return true; }
            else return false;
        }
    } ColumnCompQuality;
    struct ColumnCompColNEq_t {
        ALWAYS_INLINE bool operator() (const Query::Column& left, const Query::Column& right) {
            return 
                    (left.column ^ right.column) ||
                    (left.value ^ right.value) ||
                    (left.op ^ right.op);
            /*return  (*((uint64_t*)(&left)) ^ *((uint64_t*)&right)) ||
                    (left.value ^ right.value);*/
            //return (((uint64_t*)&left)[0]!=((uint64_t*)&right)[0]) || (((uint64_t*)&left)[1]!=((uint64_t*)&right)[1]);
            //return memcmp(&left, &right, sizeof(Query::Column));
        }
    } ColumnCompColNeq;

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
        std::vector<Query*> queries;
        ReceivedMessage *rawMsg;
        
        uint64_t validationId;
        uint64_t from,to;
        uint32_t queryCount;

        LPValidation() {}
        LPValidation(uint64_t vid, uint64_t fr, uint64_t t, uint32_t qc, ReceivedMessage* msg)
            : rawMsg(msg), validationId(vid), from(fr), to(t), queryCount(qc) {}

        ~LPValidation() {}
    }__attribute__((aligned(64)));
    struct LPValidationCompQCount {
        bool inline operator()(const LPValidation& left, const LPValidation& right){ 
            return left.queryCount > right.queryCount; 
        }
    }LPValCompQCount;

    namespace query {

        // return the new number of valid predicates
        uint32_t ALWAYS_INLINE preprocess(Query& rq) {
            if (!rq.columnCount) return 0;
            std::sort(rq.columns, rq.columns+rq.columnCount, ColumnCompCol);
            rq.columnCount = std::distance(rq.columns, 
                    std::unique(rq.columns, rq.columns+rq.columnCount, ColumnCompColEq));
            return rq.columnCount;
            //return std::distance(rq.columns, std::unique(rq.columns, rq.columns+rq.columnCount, ColumnCompColEq));
        }
        
        struct EQ{ 
            uint64_t v; uint8_t is; 
            EQ() : v(0), is(0) {} 
        };
        // return true if the query is valid otherwise false
        bool ALWAYS_INLINE preprocess(Query& rq, const size_t relCols) {
            if (rq.columnCount == 0) return true;

            //lp::simd::zero((uint8_t*)bitv, relCols*sizeof(EQ));
            EQ bitv[relCols];

            Column *qc = const_cast<Column*>(rq.columns);
            auto cb = qc, ce = cb + rq.columnCount;
            for (; cb<ce;) {
                auto& p = *cb++;
                switch (p.op) {
                    case Op::Equal:
                        if ((bitv[p.column].is) && (bitv[p.column].v != p.value)) return false;
                        bitv[p.column].is = 1; bitv[p.column].v = p.value;
                        break;
                    case Op::NotEqual:
                        if ((bitv[p.column].is) && (bitv[p.column].v == p.value)) return false;
                        break;
                    case Op::Less:
                        if ((bitv[p.column].is) && (bitv[p.column].v >= p.value)) return false;
                        break;
                    case Op::LessOrEqual:
                        if ((bitv[p.column].is) && (bitv[p.column].v > p.value)) return false;
                        break;
                    case Op::Greater:
                        if ((bitv[p.column].is) && (bitv[p.column].v <= p.value)) return false;
                        break;
                    case Op::GreaterOrEqual:
                        if ((bitv[p.column].is) && (bitv[p.column].v < p.value)) return false;
                        break;
                }
            }

            const register size_t sz = rq.columnCount;
            
            if (sz < 20) {
                if (sz == 1) return true;
                // insertion sort manual
                for (register size_t i=1; i<sz; ++i) {
                    const auto t = qc[i]; register size_t pos=i;
                    for (; pos>0; --pos) {
                        if (ColumnCompQuality(t, qc[pos-1])) {
                        //if (ColumnCompQuality(qc[pos], qc[pos-1])) {
                            //std::swap(qc[pos], qc[pos-1]);
                            qc[pos] = qc[pos-1];
                        } else break;
                    }
                    qc[pos] = t;
                }
                // unique manual
                register size_t uniq = 0;
                for (register size_t i=1; i<sz; ++i) {
                    if (ColumnCompColNeq(qc[uniq], qc[i])) {
                        ++uniq;
                        if (i > uniq) qc[uniq] = qc[i];
                        //qc[uniq] = qc[i];
                    }
                }
                rq.columnCount = uniq+1;
            } else {
                std::sort(qc, ce, ColumnCompQuality);
                rq.columnCount = std::unique(rq.columns, rq.columns+rq.columnCount, ColumnCompColEq) - qc;
            }
            
            return true;
        }
        /*
        inline uint32_t __attribute__((always_inline)) preprocess(LPQuery& lpq) {
            auto q = lpq.rawQuery;
            std::sort(q->columns, q->columns+q->columnCount, ColumnCompCol);
            auto colEnd = std::unique(q->columns, q->columns+q->columnCount, ColumnCompColEq);
            lpq.colCountUniq = std::distance(q->columns, colEnd);
            return lpq.colCountUniq;
        }
        */
        struct Satisfiability {
            uint64_t eq=UINT64_MAX, lt = UINT64_MAX, leq = UINT64_MAX, gt = 0, geq = 0;
            uint8_t EQExists;
            
            Satisfiability():eq(UINT64_MAX), lt(UINT64_MAX), leq(UINT64_MAX), gt(0), geq(0) {
                EQExists = 0;
            }
            inline void reset() {
                eq=UINT64_MAX; lt = UINT64_MAX; leq = UINT64_MAX; gt = 0; geq = 0;
                EQExists = (uint8_t)0;
            }
        };

        bool isQueryColumnUnsolvable(Column& p, Satisfiability& sat) {
            switch (p.op) {
                case Op::Equal:
                    // already found an equality check
                    if ((sat.EQExists) && (sat.eq != p.value)) return true;
                    sat.EQExists = 1;
                    sat.eq = p.value;
                    break;
                case Op::NotEqual:
                    // both equality and inequality with same value
                    if (sat.EQExists && (sat.eq == p.value)) return true;
                    break;
                case Op::Less:
                    if (sat.EQExists) {
                        if (sat.eq >= p.value) return true;
                        p.op = Op::Equal; p.value = sat.eq;
                    } else {
                        if (p.value < sat.lt) { sat.lt = p.value; sat.leq = p.value - 1; }
                    }
                    break;
                case Op::LessOrEqual:
                    if (sat.EQExists) {
                        if (sat.eq > p.value) return true;
                        p.op = Op::Equal; p.value = sat.eq;
                    } else {
                        if (p.value < sat.leq) { sat.leq = p.value; sat.lt = p.value + 1; }
                    }
                    break;
                case Op::Greater:
                    if (sat.EQExists) {
                        if (sat.eq <= p.value) return true;
                        p.op = Op::Equal; p.value = sat.eq;
                    } else {
                        if (p.value > sat.gt) { sat.gt = p.value; sat.geq = p.value + 1; }
                    }
                    break;
                case Op::GreaterOrEqual:
                    if (sat.EQExists) {
                        if (sat.eq < p.value) return true;
                        p.op = Op::Equal; p.value = sat.eq;
                    } else {
                        if (p.value > sat.geq) { sat.geq = p.value; sat.gt = p.value - 1; }
                    }
                    break;
            }

            // check for equality and constrasting ranges OR  check non-overlapping ranges
            //return false;
            return ( (sat.lt <= sat.gt) || (sat.EQExists && ((sat.eq < sat.gt) | (sat.eq > sat.lt))) );
        }

        bool isQueryUnsolvable(Column *colBegin, Column *colEnd) {
            //if (colBegin == colEnd) return false;
            Satisfiability sat;
            
            uint32_t lastCol = UINT32_MAX;
            for (; colBegin != colEnd;) {
                if (colBegin->column != lastCol) { sat.reset(); lastCol=colBegin->column; }
                if (isQueryColumnUnsolvable(*colBegin++, sat)) return true;
            }
            return false;
        }

        bool satisfiable(Query* q, uint32_t& colCountUniq) {
            if (colCountUniq == 0) return true;
            Column *qc = const_cast<Column*>(q->columns);
            auto colBegin = qc, colEnd = qc + colCountUniq; //colEnd = qc + q->columnCount;
            if (isQueryUnsolvable(colBegin, colEnd)) {
                // the query is not-satisfiable so it should be skipped-pruned
                return false;
            }
            
            //std::cerr << "\nbefore uniq: " << colCountUniq << std::endl;
            colCountUniq = std::distance(q->columns, std::unique(q->columns, q->columns+colCountUniq, ColumnCompColEq));
            colEnd = qc + colCountUniq;
            //std::cerr << "after uniq: " << colCountUniq << std::endl;
            
            //std::partial_sort(colBegin, colBegin+std::min(colCountUniq, (uint32_t)2), colEnd, ColumnCompQuality);
            std::sort(colBegin, colEnd, ColumnCompQuality);
            
            //for (auto c=colBegin;c<colEnd; ++c) std::cerr << c->column << ":" << c->op << ":" << c->value << " ";
            //std::cerr << std::endl;
            return true;
        }

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
