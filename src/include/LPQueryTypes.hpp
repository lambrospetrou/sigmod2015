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

    //static uint32_t opw[6] = { 0, 5, 1, 2, 3, 4 }; 
    //static LPOps lpopw[6] = { LPOps::Equal, LPOps::NotEqual, LPOps::Less, LPOps::LessOrEqual, LPOps::Greater, LPOps::GreaterOrEqual }; 
    //enum LPOps : uint32_t { Equal, Less, LessOrEqual, Greater, GreaterOrEqual, NotEqual }

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
            std::sort(predicates.begin(), predicates.end(), QCSortCol);
            predicates.resize(std::distance(predicates.begin(), std::unique(predicates.begin(), predicates.end(), QCEquality)));

            //if (columns.size() != columnCount) cerr << "diff: " << columnCount-columns.size() << endl;
            // reorder operators
            //for (auto& p : predicates) if (p.op == LPOps::NotEqual) p.op = LPOps::NotEqualLast ;

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
            return (left.columnCount < right.columnCount);
        }
        static bool LPQueryWeightLess(const LPQuery& left, const LPQuery& right) {
            return (left.weight < right.weight);
        }
        static bool QCEquality (const Query::Column& l, const Query::Column& r) {
            if (l.column != r.column) return false;
            if (l.op != r.op) return false;
            return l.value == r.value;
        }
    };

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

        struct PCell {
            uint64_t val;
            bool valid;
            PCell() : val(0) , valid(false) {}
        };

        struct PMatrix {
            //enum Op : uint32_t { Equal, NotEqual, Less, LessOrEqual, Greater, GreaterOrEqual };
            PCell cell[6][1000];
            PMatrix() { reset(); }
            void reset(uint32_t cols = 1000) {
                for (uint32_t op=0; op<4; ++op) {
                    for (uint32_t c=0; c<cols; ++c) {
                        cell[op][c].valid = false;
                        cell[op][c].val = UINT64_MAX;;
                    }
                }
                //memset(cell[0], 1, sizeof(PCell)*cols); // UINT64_MAX
                //memset(cell[1], 1, sizeof(PCell)*cols); // UINT64_MAX
                //memset(cell[2], 1, sizeof(PCell)*cols); // UINT64_MAX
                //memset(cell[3], 1, sizeof(PCell)*cols); // UINT64_MAX
                memset(cell[4], 0, sizeof(PCell)*cols); // 0
                memset(cell[5], 0, sizeof(PCell)*cols); // 0
            }
            PCell* operator[] (uint64_t index) { return cell[index]; }
        };

        // Return True if the query passed is fine - False if it is unsolvable
        // LPQuery& resQ : will contain the proper predicates at the end if True or will not be modified
        bool parse(const Query& q, LPQuery& resQ) {
            (void)q; (void)resQ;
            static PMatrix cell;
            cell.reset();

            PCell* EQ = cell[Op::Equal];
            PCell* LT = cell[Op::Less];
            PCell* LEQ = cell[Op::LessOrEqual];
            PCell* GT = cell[Op::Greater];
            PCell* GEQ = cell[Op::GreaterOrEqual];

            for (auto c=q.columns,cLimit=c+q.columnCount;c!=cLimit;++c) {
                const Column& p = *c;
                //std::cerr << p.column << ":" << p.op << ":" << p.value << std::endl;
                auto& pcell = cell[p.op][p.column];
                switch (p.op) {
                    case Op::Equal:
                        // already found an equality check
                        if (pcell.valid && pcell.val != p.value) return false;
                        pcell.valid = true; 
                        pcell.val = p.value; 
                        break;
                    case Op::NotEqual:
                        // both equality and inequality with same value
                        if (EQ[p.column].valid && EQ[p.column].val == p.value) return false;
                        pcell.valid = true;
                        pcell.val = p.value; // TODO - SPECIAL CARE FOR THIS - MAYBE JUST PUSH IT
                        break;
                    case Op::Less:
                        if (EQ[p.column].valid && EQ[p.column].val >= p.value) return false;
                        pcell.valid = true;
                        if (p.value < pcell.val) { pcell.val = p.value; LEQ[p.column].val = p.value - 1; }
                        break;
                    case Op::LessOrEqual:
                        if (EQ[p.column].valid && EQ[p.column].val > p.value) return false;
                        pcell.valid = true;
                        if (p.value < pcell.val) { pcell.val = p.value; LT[p.column].val = p.value + 1; }
                        break;
                    case Op::Greater:
                        if (EQ[p.column].valid && EQ[p.column].val <= p.value) return false;
                        pcell.valid = true;
                        if (p.value > pcell.val) { pcell.val = p.value; GEQ[p.column].val = p.value + 1; }
                        break;
                    case Op::GreaterOrEqual:
                        if (EQ[p.column].valid && EQ[p.column].val < p.value) return false;
                        pcell.valid = true;
                        if (p.value > pcell.val) { pcell.val = p.value; GT[p.column].val = p.value - 1; }
                        break;
                }
            }

            std::vector<Column> preds; preds.reserve(8);
            for (uint32_t op=0; op<6; ++op) {
                for (uint32_t col=0; col<1000; ++col) {
                    if (cell[op][col].valid) preds.push_back(Column(col, static_cast<Op>(op), cell[op][col].val));
                }
            }
            //for (auto& c : preds) std::cerr << c.column << ":" << c.op << ":" << c.value << std::endl;
            using std::swap;
            swap(resQ.predicates, preds);

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

        //typedef LPOps Op;
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
