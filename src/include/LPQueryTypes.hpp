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
#include <memory>


namespace lp {

    //static uint32_t opw[6] = { 0, 5, 1, 2, 3, 4 }; 
    //static LPOps lpopw[6] = { LPOps::Equal, LPOps::NotEqual, LPOps::Less, LPOps::LessOrEqual, LPOps::Greater, LPOps::GreaterOrEqual }; 
    //enum LPOps : uint32_t { Equal, Less, LessOrEqual, Greater, GreaterOrEqual, NotEqual }

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
    
    struct LPQuery {
        /// The relation
        uint32_t relationId;
        /// The number of bound columns
        uint32_t columnCount;
        /// The bindings
        std::vector<Query::Column> predicates;
        /// States whether this query can ever be true - will be set to false later by our filtering
        bool satisfiable;
        // the weight assigned to this query by the validation processor
        uint32_t weight;

        LPQuery() : relationId(-1), columnCount(0), predicates(std::vector<Query::Column>()), satisfiable(true), weight(0) {}
        LPQuery(const Query& q) : relationId(q.relationId), columnCount(q.columnCount), predicates(std::vector<Query::Column>(q.columnCount)), satisfiable(true), weight(0) {
            memcpy(predicates.data(), q.columns, sizeof(Query::Column)*columnCount);
            std::sort(predicates.begin(), predicates.end(), ColumnCompCol);
            predicates.resize(std::distance(predicates.begin(), std::unique(predicates.begin(), predicates.end(), QCEquality)));

            //if (columns.size() != columnCount) cerr << "diff: " << columnCount-columns.size() << endl;
            // reorder operators
            //for (auto& p : predicates) if (p.op == LPOps::NotEqual) p.op = LPOps::NotEqualLast ;

            columnCount = predicates.size();
        }
        // this is the default - operator< of Column
        static bool QCEquality (const Query::Column& l, const Query::Column& r) {
            if (l.column != r.column) return false;
            if (l.op != r.op) return false;
            return l.value == r.value;
        }
    };

    struct LPQueryCompSize_t {
        inline bool operator()(const LPQuery& left, const LPQuery& right) {
            return (left.columnCount < right.columnCount);
        }
    } LPQueryCompSizeLess;

    struct LPValidation {
        uint64_t validationId;
        uint64_t from,to;
        std::vector<LPQuery> queries;
        LPValidation(const ValidationQueries& v, std::vector<LPQuery> q)
            : validationId(v.validationId), from(v.from), to(v.to), queries(move(q)) {}
        LPValidation(uint64_t vid, uint64_t fr, uint64_t t, std::vector<LPQuery> q)
            : validationId(vid), from(fr), to(t), queries(q) {}
    };


    namespace query {

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

        bool parse(const Query *q, uint32_t relCols, LPQuery *nQ) {
            (void)relCols; (void)nQ;
            Column * qc = const_cast<Column*>(q->columns);
            /*
            std::cerr << std::endl;
            for (auto c=q->columns, cLimit=c+q->columnCount; c!=cLimit; ++c) {
                std::cerr << c->column << ":" << c->op << ":" << c->value << " ";
            }
            */
            //std::vector<Column> preds(q->columns, q->columns+q->columnCount);
            
            // sort the columns by column first in order to remove uniques and check satisfiability
            std::sort(qc, qc+q->columnCount, ColumnCompCol);
            auto colEnd = std::unique(qc, qc+q->columnCount, ColumnCompColEq);
            auto colBegin = qc;
            uint64_t uniqSz = std::distance(colBegin, colEnd);
            
            if (isQueryUnsolvable(colBegin, colEnd)) {
                // the query is not-satisfiable so it should be skipped-pruned
                nQ->satisfiable = false;
                nQ->columnCount = 0;
                return false;
            }
            
            // the query is satisfiable so insert it into the predicates of the passed in LPQuery
            nQ->predicates.reserve(uniqSz);
            nQ->predicates.insert(nQ->predicates.begin(), colBegin, colEnd);
            nQ->columnCount = uniqSz;
            //std::cerr << "unique: " << uend-q->columns << std::endl; 
            /*
            std::cerr << std::endl << "after unique: " << std::distance(const_cast<Column*>(q->columns), uend) << std::endl; 
            for (auto c=q->columns; c!=uend; ++c) {
            //for (auto c=uend, cLimit=const_cast<Column*>(q->columns)+q->columnCount; c!=cLimit; ++c) {
                std::cerr << c->column << ":" << c->op << ":" << c->value << " ";
            }
            if (std::distance(qc, uend) != q->columnCount) {
                std::cerr << std::endl << std::endl << " unique worked " << std::endl << std::endl; 
            }
            */
            
            //for (auto c: preds) std::cerr << c.column << ":" << c.op << ":" << c.value << " ";
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
            return false;
        }
    }
}

#endif
