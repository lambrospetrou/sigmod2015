// SIGMOD Programming Contest 2015 - reference implementation
//
// This code is intended to illustrate the given challenge, it
// is not particular efficient.
//---------------------------------------------------------------------------
// This is free and unencumbered software released into the public domain.
//
// Anyone is free to copy, modify, publish, use, compile, sell, or
// distribute this software, either in source code form or as a compiled
// binary, for any purpose, commercial or non-commercial, and by any
// means.
//
// In jurisdictions that recognize copyright laws, the author or authors
// of this software dedicate any and all copyright interest in the
// software to the public domain. We make this dedication for the benefit
// of the public at large and to the detriment of our heirs and
// successors. We intend this dedication to be an overt act of
// relinquishment in perpetuity of all present and future rights to this
// software under copyright law.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
//
// For more information, please refer to <http://unlicense.org/>
//---------------------------------------------------------------------------

#include "include/CIndex.hpp"
#include "include/DoubleIter.hpp"
#include "include/aligned_allocator.hpp"
#include "include/LPUtils.hpp"
#include "include/ReaderIO.hpp"
#include "include/LPThreadpool.hpp"
#include "include/SingleTaskPool.hpp"

#include <iostream>
#include <fstream>
#include <ios>
#include <cstdio>
#include <iterator>
#include <map>
#include <unordered_map>
#include <set>
#include <vector>
#include <bitset>
#include <cstdint>
#include <algorithm>
#include <utility>
#include <functional>
#include <cstring>
#include <cstdlib>
#include <climits>
#include <sys/time.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <memory>
#include <system_error>
#include <exception>
#include <functional>

#include "include/atomicwrapper.hpp"
#include "include/LPTimer.hpp"
#include "include/BoundedQueue.hpp"
#include "include/LPSpinLock.hpp"

#include "include/cpp_btree/btree_map.h"

#include "include/ReferenceTypes.hpp"
#include "include/LPQueryTypes.hpp"


//---------------------------------------------------------------------------
using namespace std;
using namespace lp;

#define LPDEBUG  

//---------------------------------------------------------------------------
// Begin reference implementation
//---------------------------------------------------------------------------

template<typename T>
ostream& operator<< (ostream& os, const vector<T> v) {
    if (v.empty()) return os;
    os << *v.begin();
    for (auto it=v.begin()+1,vend=v.end(); it!=vend; ++it) {
        os << " " << *it;
    }
    return os;
}
std::ostream& operator<< (std::ostream& os, const Query::Column& o) {
    os << "[" << o.column << ":" << o.op << ":" << o.value << "]";
    return os;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

typedef uint64_t* tuple_t;
typedef Query::Column::Op  Op;

#define CACHE_LINE_SIZE 64
#define CACHE_ALIGNMENT 64

#define ALIGNED_DATA __attribute__((aligned(CACHE_ALIGNMENT)))

template<typename T>
//using vector_a = std::vector<T, aligned_allocator<T, 16>>;
using vector_a = std::vector<T>;

// Custom data structures to hold data
struct CTransStruct {
    uint64_t value;
    tuple_t tuple;
    CTransStruct (uint64_t v, tuple_t t) : value(v), tuple(t) {}
};
struct CTRSValueLessThan_t {
    ALWAYS_INLINE bool operator()(const CTransStruct& l, const CTransStruct& r) {
        return l.value < r.value;
    }
    ALWAYS_INLINE bool operator() (const CTransStruct& o, uint64_t target) {
        return o.value < target;
    }
    ALWAYS_INLINE bool operator() (uint64_t target, const CTransStruct& o) {
        return target < o.value;
    }
} ColTransValueLess;
std::ostream& operator<< (std::ostream& os, const CTransStruct& o) {
    os << "{" << "-" << o.value << "}";
    return os;
}


// transactions in each relation column - all tuples of same transaction in one vector
struct Metadata_t {
    uint32_t tpl_id;
    tuple_t tuple;
};
struct ColumnTransaction_t {
    vector_a<uint64_t> values;
    vector_a<Metadata_t> tuples;
    uint64_t trans_id;
    ColumnTransaction_t(uint64_t tid) : values(vector_a<uint64_t>()), tuples(vector_a<Metadata_t>()), trans_id(tid) {}
} ALIGNED_DATA;
struct ColumnStruct {
    // the trans_id the transactions are updated to inclusive
    //vector<ColumnTransaction_t> transactions;
    //vector<uint64_t> transactionsORs;
    CIndex transactions;
    uint64_t transTo;

    //char padding[8]; //for false sharing

    ColumnStruct() : transTo(0) {}
} ALIGNED_DATA;
struct CTRSLessThan_t {
    ALWAYS_INLINE bool operator() (const ColumnTransaction_t& left, const ColumnTransaction_t& right) {
        return left.trans_id < right.trans_id;
    }
    ALWAYS_INLINE bool operator() (const ColumnTransaction_t& o, uint64_t target) {
        return o.trans_id < target;
    }
    ALWAYS_INLINE bool operator() (uint64_t target, const ColumnTransaction_t& o) {
        return target < o.trans_id;
    }
} CTRSLessThan;


struct RelationColumns {
    //ColumnStruct *columns;
    std::unique_ptr<ColumnStruct[]> columns;
};
static std::unique_ptr<RelationColumns[]> gRelColumns;
static vector<uint32_t> gRequiredColumns;

//////////////////////////////////////////////////////////

// QUERIES INDEX

struct QMeta_t{
    uint64_t from;
    uint64_t to;
    Query *rq;
    LPValidation *lpv;
    uint64_t value;
};
struct QMVLess_t {
    ALWAYS_INLINE bool operator() (const QMeta_t& left, const QMeta_t& right) {
        return left.value < right.value;
    }
    ALWAYS_INLINE bool operator() (const QMeta_t& o, uint64_t target) {
        return o.value < target;
    }
    ALWAYS_INLINE bool operator() (uint64_t target, const QMeta_t& o) {
        return target < o.value;
    }
} QMVLess;
struct CQ_t {
    vector<QMeta_t> queries;
    std::mutex mtx;
    //LPSpinLock mtx;
} ALIGNED_DATA;
struct RelQ_t {
    std::unique_ptr<CQ_t[]> columns;
};
static std::unique_ptr<RelQ_t[]> gRelQ;

static std::atomic<uint64_t> gNextEQCol;
static vector<uint64_t> gEQCols;

//////////////////////////////////////////////////////////
struct RelTransLog {
    uint64_t trans_id;
    uint64_t last_del_id;
    uint64_t aliveTuples;
    uint64_t rowCount;
    //uint64_t *tuples;
    std::unique_ptr<uint64_t[]> tuples;

    char padding[22]; //for false sharing

    RelTransLog(uint64_t tid, uint64_t* tpl, uint64_t _rowCount) : trans_id(tid), last_del_id(tid), aliveTuples(_rowCount), rowCount(_rowCount) {
        if (tpl != nullptr) tuples.reset(tpl);
    } 
} ALIGNED_DATA;

struct RTLComp_t {
    inline bool operator() (const RelTransLog& l, const RelTransLog& r) {
        return l.trans_id < r.trans_id;
    }
    inline bool operator() (const RelTransLog& o, uint64_t target) {
        return o.trans_id < target;
    }
    inline bool operator() (const unique_ptr<RelTransLog>& l, const unique_ptr<RelTransLog>& r) {
        return l->trans_id < r->trans_id;
    }
    inline bool operator() (const unique_ptr<RelTransLog>& o, uint64_t target) {
        return o->trans_id < target;
    }
} RTLComp;

struct RelationStruct {
    vector<pair<uint64_t, vector<tuple_t>>> transLogTuples;
    vector<unique_ptr<RelTransLog>> transLog;
    btree::btree_map<uint32_t, pair<uint64_t, uint64_t*>> insertedRows;

    //CIndex primaryIndex;

    char padding[8]; //for false sharing
} ALIGNED_DATA;
struct TransLogComp_t {
    inline bool operator()(const pair<uint64_t, vector<tuple_t>>& l, const uint64_t target) {
        return l.first < target;
    }
    inline bool operator()(const uint64_t target, const pair<uint64_t, vector<tuple_t>>& r) {
        return target < r.first;
    }
} TransLogComp;

// general relations
static std::unique_ptr<RelationStruct[]> gRelations;

struct TRMapPhase {
    uint64_t trans_id;
    uint64_t *values; // delete op => row keys to delete | insert => tuples
    uint32_t rowCount;
    bool isDelOp;

    TRMapPhase(uint64_t tid, bool isdel, uint32_t rows, uint64_t *vals)
        : trans_id(tid), values(vals), rowCount(rows), isDelOp(isdel) {
        }
} ALIGNED_DATA;
struct TransMapPhase_t {
    inline bool operator()(const TRMapPhase& l, const TRMapPhase& r) {
        if (l.trans_id < r.trans_id) return true;
        else if (r.trans_id < l.trans_id) return false;
        else return l.isDelOp && !r.isDelOp;
    }
} TRMapPhaseByTrans;

typedef std::mutex RelTransLock;
//typedef LPSpinLock RelTransLock;
static std::unique_ptr<RelTransLock[]> gRelTransMutex;
static std::unique_ptr<vector<TRMapPhase>[]> gTransParseMapPhase;

//////////////////////////////////////////////////////////////

static uint32_t NUM_RELATIONS;
static std::unique_ptr<uint32_t[]> gSchema;

///////// AUXILIARY STRUCTURES FOR THE WHOLE PROGRAM

static vector<LPValidation> gPendingValidations;
static std::mutex gPendingValidationsMutex;
//static LPSpinLock gPendingValidationsMutex;

//typedef atomic_wrapper<bool> PendingResultType;
typedef char PendingResultType;
static vector<PendingResultType> gPendingResults;
static vector<pair<uint64_t,bool>> gQueryResults;
static uint64_t gPVunique;


////////////////////////////////////////////////////

LPTimer_t LPTimer;

struct GlobalState {
    enum State : uint32_t { SCHEMA, TRANSACTION, VALIDATION, FORGET, FLUSH };
    State state;
    uint32_t nThreads;
} Globals;

//---------------------------------------------------------------------------

// JUST SOME FUNCTION DECLARATIONS THAT ARE DEFINED BELOW
static void checkPendingValidations(ISingleTaskPool*, ISingleTaskPool*);
void processPendingValidationsTask(uint32_t nThreads, uint32_t tid, void *args);

static void processTransactionMessage(const Transaction& t, ReceivedMessage *msg); 
void processPendingIndexTask(uint32_t nThreads, uint32_t tid, void *args);

void processUpdateIndexTask(uint32_t nThreads, uint32_t tid, void *args);
static void updateRelCol(uint32_t tid, uint32_t ri, uint32_t col);
///--------------------------------------------------------------------------
///--------------------------------------------------------------------------

static void processDefineSchema(const DefineSchema& d) {
    gSchema.reset(new uint32_t[d.relationCount]);
    memcpy(gSchema.get(), d.columnCounts, sizeof(uint32_t)*d.relationCount);
    NUM_RELATIONS = d.relationCount;
    //cerr << "relations: " << NUM_RELATIONS << endl;

    gRelTransMutex.reset(new RelTransLock[d.relationCount]);
    gTransParseMapPhase.reset(new vector<TRMapPhase>[d.relationCount]);

    gRelations.reset(new RelationStruct[d.relationCount]);
    gRelColumns.reset(new RelationColumns[d.relationCount]);
    gRelQ.reset(new RelQ_t[d.relationCount]);
    //cerr << endl << "relations: " << NUM_RELATIONS << endl;
    const uint32_t rels = d.relationCount;
    for(uint32_t ri=0; ri<rels; ++ri) {
        //cerr << " " << gSchema[ci];
        gRelColumns[ri].columns.reset(new ColumnStruct[gSchema[ri]]);
        gRelQ[ri].columns.reset(new CQ_t[gSchema[ri]]);
        const uint32_t colsz = gSchema[ri];
        for (uint32_t ci=0; ci<colsz; ++ci) {
            gRequiredColumns.push_back(lp::validation::packRelCol(ri, ci));
            gEQCols.push_back(lp::validation::packRelCol(ri, ci));
        }
    }
}
//---------------------------------------------------------------------------

#ifdef LPDEBUG
static uint64_t gTotalTransactions = 0, gTotalTuples = 0, gTotalValidations = 0;
#endif
//---------------------------------------------------------------------------
static void ALWAYS_INLINE processValidationQueries(const ValidationQueries& v, ReceivedMessage *msg, uint64_t tid = 0) {
#ifdef LPDEBUG
    auto start = LPTimer.getChrono();    
#endif 
    (void)tid;

    //cerr << v.validationId << "====" << v.from << ":" << v.to << "=" << v.queryCount << endl;
    gPendingValidations.emplace_back(v.validationId, v.from, v.to, v.queryCount, msg);    
    // update the global pending validations to reflect this new one
    ++gPVunique;

#ifdef LPDEBUG
    LPTimer.validations += LPTimer.getChrono(start);
#endif
    return;
}

//---------------------------------------------------------------------------
static void processFlush(const Flush& f, bool isTestdriver) {
#ifdef LPDEBUG
    auto start = LPTimer.getChrono();
#endif
    static char zero = '0';
    static char one = '1';
    if (!gQueryResults.empty()) {
        uint64_t removed = 0;
        for (auto& vp : gQueryResults) {
            if (vp.first > f.validationId) break;
            if (unlikely(!isTestdriver))
                cout << (vp.second == true ? one : zero); 
            ++removed;  
        }
        if (removed > 0) {
            if (removed == gQueryResults.size()) gQueryResults.clear();
            else gQueryResults.erase(gQueryResults.begin(), gQueryResults.begin()+removed);
        }
        cout.flush();
    }
#ifdef LPDEBUG
    LPTimer.flushes += LPTimer.getChrono(start);
#endif
}
//---------------------------------------------------------------------------


static void ALWAYS_INLINE forgetRel(uint64_t trans_id, uint32_t ri) {
    auto& cRelCol = gRelColumns[ri];
    auto& transLogTuples = gRelations[ri].transLogTuples;

    if (transLogTuples.empty() || transLogTuples[0].first > trans_id) return;

    // delete the transLogTuples
    transLogTuples.erase(transLogTuples.begin(), 
            upper_bound(transLogTuples.begin(), transLogTuples.end(), trans_id,
                [](const uint64_t target, const pair<uint64_t, vector<tuple_t>>& o){ return target < o.first; })
            );

    //cerr << ri << " = " << cRelCol.columns[0].transactions.size() << " | " << cRelCol.columns[1].transactions.size() << endl;
    for (uint32_t ci=0,sz=gSchema[ri]; ci<sz; ++ci) {
        cRelCol.columns[ci].transactions.erase(trans_id);
    }
    // clean the transactions log
}

static atomic<uint64_t> gNextFRel;
void processForgetThreaded(uint32_t nThreads, uint32_t tid, void *args) { 
    (void)tid; (void)nThreads; (void)args;// to avoid unused warning
    auto f = *reinterpret_cast<Forget*>(args);
    for (uint64_t ri = gNextFRel++; ri < NUM_RELATIONS; ri=gNextFRel++) { 
        forgetRel(f.transactionId, ri);
    }
}

static void processForget(const Forget& f, ISingleTaskPool* pool) {(void)pool;
#ifdef LPDEBUG
    auto start = LPTimer.getChrono();
#endif
    /*
       gNextFRel = 0; 
       pool->startSingleAll(processForgetThreaded, (void*)&f);
       pool->waitSingleAll();
     */
    // delete the transactions from the columns index
    for (uint32_t i=0; i<NUM_RELATIONS; ++i) {
        forgetRel(f.transactionId, i);
    }
#ifdef LPDEBUG
    LPTimer.forgets += LPTimer.getChrono(start);
#endif
}
//---------------------------------------------------------------------------



/**
  * THE FOLLOWING STRUCTURE IS USED IN ORDER TO ALLOW OVERLAPPING AMONG 
  * THE EXECUTION PHASES. 
  * TODO refactor ? - Last minute, so just leave it here,
  */ 
struct FR_t{
    atomic_wrapper<bool> f1;
    atomic_wrapper<bool> f2;
    atomic_wrapper<uint64_t> colsdone;
    atomic_wrapper<uint64_t> cols;
    FR_t() : f1(false), f2(false), colsdone(0), cols(0) {}
    void reset() {
        f1 = false; 
        f2 = false; 
        colsdone = 0; 
        cols = 0;
    }
};
static vector<FR_t> gFR;
static std::mutex mMtxF1;
static std::condition_variable mCondF1;
static std::atomic<uint32_t> finishedF1;

void resetUpdateStats() {
    for (uint32_t ri=0; ri<NUM_RELATIONS; ri++) {
        gFR[ri].reset();
    }
    finishedF1 = 0;
}

//---------------------------------------------------------------------------
/////////////////// MAIN-READING STRUCTURES ///////////////////////

int main(int argc, char**argv) {
    uint64_t numOfThreads = 1;
    if (argc > 1) {
        numOfThreads = strtol(argv[1], NULL, 10);
        if (numOfThreads == LONG_MAX)
            numOfThreads = 1;
    }
    //cerr << "Number of threads: " << numOfThreads << endl;
    Globals.nThreads = numOfThreads;

    std::ifstream ifs; bool isTestdriver = false;  
    ReaderIO* msgReader;
    if (argc > 2) {
        // there is a file argument os use the test driver reader
        ifs.open(argv[2], ios::in | ios::binary);
        if (!ifs) {
            cerr << "ERROR opening the file passed as 2nd argument!!!" << endl;
            abort();
        }
        isTestdriver = true;
        //msgReader = ReaderIOFactory::createAsync(ifs, true);
        msgReader = ReaderIOFactory::create(ifs, true);
    } else { 
        msgReader = ReaderIOFactory::createAsync(stdin);
        //msgReader = ReaderIOFactory::create(stdin);
        //msgReader = ReaderIOFactory::create(cin);
    }

    // do some initial reserves or initializations
    gPendingValidations.reserve(512); 
    for (uint32_t i=0; i<NUM_RELATIONS; ++i) {
        gTransParseMapPhase[i].reserve(512);
        gRelations[i].transLog.reserve(512);
        gRelations[i].transLogTuples.reserve(512);
    }

    // allocate the workers
    SingleTaskPool workerThreads(numOfThreads, processPendingValidationsTask);
    //SingleTaskPool workerThreads(2, processPendingValidationsTask);
    workerThreads.initThreads();

    cerr << "ColumnStruct: " << sizeof(ColumnStruct) << " RelTransLog: " << sizeof(RelTransLog) << " RelationStruct: " << sizeof(RelationStruct) << " CTransStruct: " << sizeof(CTransStruct) << endl;

    try {

#ifdef LPDEBUG
        uint64_t totalForgets = 0, totalFlushes = 0;
#endif

        while (true) {
#ifdef LPDEBUG 
            auto startInner = LPTimer.getChrono();
#endif
            ReceivedMessage *msg = msgReader->nextMsg();
#ifdef LPDEBUG
            LPTimer.reading += LPTimer.getChrono(startInner);
#endif 

            auto& head = msg->head;
            switch (head.type) {
                case MessageHead::ValidationQueries: 
                    {    Globals.state = GlobalState::VALIDATION;
#ifdef LPDEBUG
                        ++gTotalValidations; // this is just to count the total validations....not really needed!
#endif
                        processValidationQueries(*reinterpret_cast<const ValidationQueries*>(msg->data.data()), msg); 

                        break;
                    }
                case MessageHead::Transaction: 
#ifdef LPDEBUG
                    ++gTotalTransactions; 
#endif
                    Globals.state = GlobalState::TRANSACTION;
                    processTransactionMessage(*reinterpret_cast<const Transaction*>(msg->data.data()), msg); 
                    break;
                case MessageHead::Flush:  
#ifdef LPDEBUG
                    ++totalFlushes; 
#endif
                    // check if we have pending transactions to be processed
                    checkPendingValidations(&workerThreads, nullptr);
                    Globals.state = GlobalState::FLUSH;
                    processFlush(*reinterpret_cast<const Flush*>(msg->data.data()), isTestdriver); 
                    delete msg;
                    break;
                case MessageHead::Forget: 
#ifdef LPDEBUG
                    ++totalForgets; 
#endif
                    // check if we have pending transactions to be processed
                    checkPendingValidations(&workerThreads, nullptr);
                    Globals.state = GlobalState::FORGET;
                    processForget(*reinterpret_cast<const Forget*>(msg->data.data()), &workerThreads); 
                    delete msg;
                    break;
                case MessageHead::DefineSchema: 
                    Globals.state = GlobalState::SCHEMA;
                    processDefineSchema(*reinterpret_cast<const DefineSchema*>(msg->data.data()));
                    delete msg;
                    break;
                case MessageHead::Done: 
                    {
#ifdef LPDEBUG
                        cerr << "  :::: " << LPTimer << endl << "total validations: " << gTotalValidations << " trans: " << gTotalTransactions << " tuples: " << gTotalTuples << " forgets: " << totalForgets << " flushes: " << totalFlushes <<  endl; 
#endif              
                        workerThreads.destroy();
                        delete msgReader;
                        return 0;
                    }
                default: cerr << "malformed message" << endl; abort(); // crude error handling, should never happen
            }

        }
    } catch (const std::exception& e) { cerr <<  "exception " <<  e.what() << endl; }
}
//---------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////////////////
/*
   struct TRMapPhase {
   uint64_t trans_id;
   bool isDelOp;
   uint32_t rowCount;
   uint64_t *values; // delete op => row keys to delete | insert => tuples
   }
//static unique_ptr<vector<TRMapPhase>[]> gTransParseMapPhase;
 */
static void processTransactionMessage(const Transaction& t, ReceivedMessage *msg) {
#ifdef LPDEBUG
    auto start = LPTimer.getChrono();
#endif
    const char* reader=t.operations;
    // Delete all indicated tuples
    for (uint32_t index=0;index!=t.deleteCount;++index) {
        auto& o=*reinterpret_cast<const TransactionOperationDelete*>(reader);
        {// start of lock_guard
            uint64_t *ptr = new uint64_t[o.rowCount];
            for (uint32_t c=0; c<o.rowCount; ++c) ptr[c] = o.keys[c]; // faster sometimes
            //memcpy(ptr, o.keys, sizeof(uint64_t)*o.rowCount);
            gTransParseMapPhase[o.relationId].emplace_back(t.transactionId, true, o.rowCount, ptr);
#ifdef LPDEBUG
            gTotalTuples += o.rowCount; // contains more than the true
#endif
        }// end of lock_guard
        // advance to the next Relation deletions
        reader+=sizeof(TransactionOperationDelete)+(sizeof(uint64_t)*o.rowCount);
    }

    // Insert new tuples
    for (uint32_t index=0;index!=t.insertCount;++index) {
        auto& o=*reinterpret_cast<const TransactionOperationInsert*>(reader);
        const uint32_t relCols = gSchema[o.relationId];
        {// start of lock_guard
            uint64_t* tptr = new uint64_t[relCols*o.rowCount];
            uint32_t sz = relCols*o.rowCount;
            for (uint32_t c=0; c<sz; ++c) tptr[c] = o.values[c];
            //memcpy(tptr, o.values, sizeof(uint64_t)*relCols*o.rowCount);
            gTransParseMapPhase[o.relationId].emplace_back(t.transactionId, false, o.rowCount, tptr);
#ifdef LPDEBUG
            ++gTotalTuples;
#endif
        }// end of lock_guard
        // advance to next Relation insertions
        reader+=sizeof(TransactionOperationInsert)+(sizeof(uint64_t)*o.rowCount*relCols);
    }

    delete msg;
#ifdef LPDEBUG
    LPTimer.transactions += LPTimer.getChrono(start);
#endif
}


static std::atomic<uint64_t> gNextIndex;  // for the index
static std::atomic<uint32_t> gNextReqCol; // for the update columns

void processPendingIndexTask(uint32_t nThreads, uint32_t tid, void *args) {
    (void)tid; (void)nThreads; (void)args;// to avoid unused warning
    //cerr << "::: tid " << tid << "new" << endl;
    for (uint64_t ri = gNextIndex++; ri < NUM_RELATIONS; ri=gNextIndex++) {

        // take the vector with the transactions and sort it by transaction id in order to apply them in order
        auto& relTrans = gTransParseMapPhase[ri];
        if (unlikely(relTrans.empty())) { 
            gFR[ri].f1 = true; gFR[ri].f2 = true;
            ++finishedF1; 
            continue; 
        }

        auto& relation = gRelations[ri];
        uint32_t relCols = gSchema[ri];

        uint64_t lastTransId = relTrans[0].trans_id;

        // for each transaction regarding this relation
        vector<tuple_t> operations;
        operations.reserve(64);
        for (auto& trans : relTrans) {
            if (trans.trans_id != lastTransId) {
                // store the tuples for the last transaction just finished
                if (likely(!operations.empty()))
                    relation.transLogTuples.emplace_back(lastTransId, operations);
                lastTransId = trans.trans_id;
                operations.resize(0);
            }

            if (trans.isDelOp) {
                // this is a delete operation
                for (const uint64_t* key=trans.values,*keyLimit=key+trans.rowCount;key!=keyLimit;++key) {
                    auto lb = relation.insertedRows.find(*key);
                    if (lb != relation.insertedRows.end()) {
                        // update the relation transactions - transfer ownership of the tuple
                        operations.push_back(lb->second.second);
                        // remove the row from the relations table 
                        relation.insertedRows.erase(lb);
                    }
                }
                delete[] trans.values;
            } else {
                // this is an insert operation
                for (const uint64_t* values=trans.values,*valuesLimit=values+(trans.rowCount*relCols);values!=valuesLimit;values+=relCols) {
                    tuple_t vals = const_cast<uint64_t*>(values);
                    operations.push_back(vals);
                    // finally add the new tuple to the inserted rows of the relation
                    relation.insertedRows[values[0]]=move(make_pair(trans.trans_id, vals));
                }
                // TODO - THIS HAS TO BE IN ORDER - each relation will have its own transaction history from now on
                relation.transLog.emplace_back(new RelTransLog(trans.trans_id, trans.values, trans.rowCount));
            }
        }
        // store the last transaction data
        // store the operations for the last transaction
        if (likely(!operations.empty()))
            relation.transLogTuples.emplace_back(lastTransId, move(operations));

        // F1 - Phase 1 - finished for relation 'ri'
        gFR[ri].f1 = true;
        ++finishedF1; 
    } // end of all relations
}

static void updateRelCol(uint32_t tid, uint32_t ri, uint32_t col) { (void)tid;
    auto& relation = gRelations[ri];
    auto& relColumn = gRelColumns[ri].columns[col];

    uint64_t updatedUntil = relColumn.transTo;
    if (relation.transLogTuples.empty() || (updatedUntil > relation.transLogTuples.back().first)) return;

    // Use lower_bound to automatically jump to the transaction to start
    //auto transFrom = lower_bound(relation.transLogTuples.begin(), relation.transLogTuples.end(), updatedUntil, TransLogComp);
    //auto tEnd=relation.transLogTuples.end();
    auto tEnd=relation.transLogTuples.data()+relation.transLogTuples.size();
    auto transFrom = lower_bound(relation.transLogTuples.data(), tEnd, updatedUntil, TransLogComp);

    auto& cindex = relColumn.transactions;
    // for all the transactions in the relation
    for(auto trp=transFrom; trp<tEnd; ++trp) {
        // allocate vectors for the current new transaction to put its data
        auto trans_id = trp->first;
        CIndex::Bucket& trb = *cindex.bucketNext(trans_id, col==0);
       
        const size_t trpsz = trp->second.size();
        auto ptrs = trb.resizeAndGetPtr(trpsz);
        auto tplPtr = ptrs.second;
        for (auto tpl : trp->second) { 
             *tplPtr++ = {tpl[col], trans_id, tpl}; 
        }
    }
    // no need to check for empty since now we update all the columns and there is a check for emptyness above
    relColumn.transTo = relation.transLogTuples.back().first + 1;
    cindex.sortFrom(transFrom->first, col == 0);
    //if (col == 0) cindex.sortFrom(transFrom->first, true);
    //else cindex.sortFrom(transFrom->first);
}

/**
  * QUERY INVERTED INDEX
  */
static uint64_t resIndexOffset = 0;
static std::atomic<uint64_t> gNextPending;
void createQueryIndexTask(uint32_t nThreads, uint32_t tid, void *args) { (void)nThreads; (void)tid; (void)args;
    uint64_t totalPending = gPendingValidations.size();
    // get a validation ID - atomic operation
    for (uint64_t vi = gNextPending++; vi < totalPending; vi=gNextPending++) {
        auto& v = gPendingValidations[vi];
        const ValidationQueries& vq = *reinterpret_cast<ValidationQueries*>(v.rawMsg->data.data());

        //cerr << "\n\t----- VAL: " << vq.validationId << " -----" << endl;
        const char* qreader = vq.queries; uint32_t columnCount;
        for (uint32_t i=0; i<vq.queryCount; ++i, qreader+=sizeof(Query)+(sizeof(Query::Column)*columnCount)) {
            Query *rq=const_cast<Query*>(reinterpret_cast<const Query*>(qreader));
            columnCount = rq->columnCount;

            if (columnCount == 0) { v.queries.push_back(rq); continue; }

            if (!lp::query::preprocess(*rq, gSchema[rq->relationId])) { continue; }
            auto pFirst = (Query::Column)rq->columns[0];
            if (pFirst.op == Op::Equal) {
                /*
                   struct QMeta_t{
                   uint64_t from;
                   uint64_t to;
                   Query *rq;
                   LPValidation *lpv;
                   uint64_t value;
                   };
                 */
                auto& relcol = gRelQ[rq->relationId].columns[pFirst.column];
                relcol.mtx.lock();
                relcol.queries.push_back({vq.from, vq.to, rq, &v, pFirst.value});
                relcol.mtx.unlock();
            } else {
                v.queries.push_back(rq);
            }
            //v.queries.push_back(rq);
        }
    } // end for all validations
}

// update the column-index with busy waiting at the end
void processUpdateIndexTask(uint32_t nThreads, uint32_t tid, void *args) {
    (void)tid; (void)nThreads; (void)args;// to avoid unused warning
    for (size_t ri=0; ri<NUM_RELATIONS; ++ri) {
        if (!gFR[ri].f1) continue;
        auto& rctodo = gFR[ri].cols;
        size_t relCols = gSchema[ri];
        for (uint64_t rc = rctodo.get()++; rc < relCols; rc=rctodo.get()++) {
            updateRelCol(tid, ri, rc);
            // we are the last thread to finish
            if (++gFR[ri].colsdone.get() == relCols) {
                gFR[ri].f2 = true;
            }
        }
    } // end of while columns to update     
    
    bool allfinished = false;
    while (!allfinished) {
        for (size_t ri=0; ri<NUM_RELATIONS; ++ri) {
            if (!gFR[ri].f1) { allfinished = false; continue; }
            if (gFR[ri].f2) { continue; } // already finished

            auto& rctodo = gFR[ri].cols;
            size_t relCols = gSchema[ri];
            for (uint64_t rc = rctodo.get()++; rc < relCols; rc=rctodo.get()++) {
                updateRelCol(tid, ri, rc);
                // we are the last thread to finish
                if (++gFR[ri].colsdone.get() == relCols) {
                    gFR[ri].f2 = true;
                }
            }
        }
        
        if (finishedF1 < NUM_RELATIONS) lp_spin_sleep(std::chrono::microseconds(10));
        else { allfinished = true; }
    } // all relations finished
}

void parallelTask1(uint32_t nThreads, uint32_t tid, void *args) { (void)nThreads; (void)tid; (void)args;
    processPendingIndexTask(nThreads, tid, args); // F1
    createQueryIndexTask(nThreads, tid, args);    // F3
    processUpdateIndexTask(nThreads, tid, args);  // F2
}

///////////////////////////////////////////////////////
// THE MAIN PROCESSING METHOD - CALLS ALL 5 STEPS
static void checkPendingValidations(ISingleTaskPool *pool, ISingleTaskPool *pool2) { (void)pool2;
    if (unlikely(gPendingValidations.empty())) return;
    if (unlikely(gFR.empty())) gFR.resize(NUM_RELATIONS);
#ifdef LPDEBUG
    auto startQuery = LPTimer.getChrono();
#endif
    
    resetUpdateStats();
    // trans-index & qindex
    gNextPending = 0; // F3
    gNextIndex = 0; // F1
    gNextReqCol = 0; // F2 
    pool->startSingleAll(parallelTask1); // F1 + F3
    pool->waitSingleAll(); 

    // this is needed by the trans-index after it has finished - F1
    for (uint32_t r=0; r<NUM_RELATIONS; ++r) gTransParseMapPhase[r].clear();

#ifdef LPDEBUG
    LPTimer.queryIndex += LPTimer.getChrono(startQuery);
#endif

////////////////////////////// ------------ VALIDATIONS ------------- //////////////////
#ifdef LPDEBUG
    auto start = LPTimer.getChrono();
#endif

    // find the MIN validation ID to coordinate the indexing of the results
    resIndexOffset = gPendingValidations[0].validationId; // only master handles the messages now so they are in order

    const size_t gPRsz = gPendingResults.size();
    if (gPVunique > gPRsz)
        gPendingResults.resize(gPVunique);
    for (auto gpr=gPendingResults.data(), gpre=gpr+gPRsz; gpr<gpre; ) *gpr++ = 0;

    gNextPending = 0; gNextEQCol = 0;
    pool->startSingleAll(processPendingValidationsTask);
    pool->waitSingleAll();

    // update the results - you can get the validation id by adding resIndexOffset to the position
    for (uint64_t i=0, valId=resIndexOffset; i<gPVunique; ++i, ++valId) { 
        gQueryResults.emplace_back(valId, gPendingResults[i]);
    }
    for (auto& v : gPendingValidations) delete v.rawMsg;
    gPendingValidations.clear();
    gPVunique = 0;
    
    // clear the inverted query index 
    uint32_t ri, ci;
    for (uint64_t rc = 0, totalCols = gEQCols.size(); rc < totalCols; ++rc) {
        lp::validation::unpackRelCol(gEQCols[rc], ri, ci);
        gRelQ[ri].columns[ci].queries.clear();
    }

#ifdef LPDEBUG
    LPTimer.validationsProcessing += LPTimer.getChrono(start);
#endif
}

//////////////////////////////////////////////////////////
// VALIDATION
//////////////////////////////////////////////////////////

typedef Metadata_t TupleType;
typedef Query::Column* PredIter;

bool ALWAYS_INLINE isTupleConflict(PredIter cbegin, PredIter cend, const tuple_t& tup) {
    for (; cbegin<cend;) {
        register auto& c = *cbegin++;
        // make the actual check
        switch (c.op) {
            case Op::Equal: 
                if((tup[c.column] != c.value)) return false; 
                break;
            case Op::Less: 
                if((tup[c.column]>=c.value)) return false; 
                break;
            case Op::LessOrEqual: 
                if((tup[c.column]>c.value)) return false; 
                break;
            case Op::Greater: 
                if((tup[c.column]<=c.value)) return false; 
                break;
            case Op::GreaterOrEqual: 
                if((tup[c.column]<c.value)) return false; 
                break;
            case Op::NotEqual: 
                if(tup[c.column] == c.value) return false; 
                break;
        } 
    } // end of single query predicates
    return true;    
}

template<typename Tuple_t>
bool ALWAYS_INLINE isTupleRangeConflict(Tuple_t *tupFrom, Tuple_t *tupTo, PredIter cbegin, PredIter cend) {
    if (cbegin == cend && tupTo-tupFrom != 0) return true;
    for(; tupFrom!=tupTo; ++tupFrom) {
        if (isTupleConflict(cbegin, cend, tupFrom->tuple)) return true;
    } // end of all tuples for this transaction
    return false;
    //return std::find_if(tupFrom, tupTo, [=](TupleType& tup) { return isTupleConflict(cbegin, cend, tup);}) != tupTo;
}

template<typename Tuple_t>
bool ALWAYS_INLINE isTupleRangeConflict(Tuple_t *tupFrom, Tuple_t *tupTo, PredIter cbegin, PredIter cend, uint64_t vfrom, uint64_t vto) {
    const uint64_t rangediff = vto - vfrom;
    for(; tupFrom!=tupTo; ++tupFrom) {
        if (((uint64_t)(tupFrom->trans_id - vfrom) <= (rangediff)) 
                && isTupleConflict(cbegin, cend, tupFrom->tuple)) return true;
    } // end of all tuples for this transaction
    return false;
}
    
/*

// THESE FUNCTIONS IN COMMENT WERE USED BEFORE THE INDEX WAS CONVERTED TO BUCKETS
// 2 days before the end of the contest.

static bool inline isTransactionConflict(const ColumnTransaction_t& transaction, Column pFirst) {
    //cerr << pFirst << " sz: " << transValues.size() << " " << transValues[0].value << ":" << transValues.back().value <<  endl;
    auto& transValues = transaction.values;
    switch (pFirst.op) {
        case Op::Equal: 
            //return lp::utils::binary_cmov(transValues.data(), transValues.size(), pFirst.value); 
            return std::binary_search(transValues.data(), transValues.data() + transValues.size(), pFirst.value); 
        case Op::Less: 
            return transValues[0] < pFirst.value;                   
        case Op::LessOrEqual: 
            return transValues[0] <= pFirst.value;                   
        case Op::Greater: 
            return transValues.back() > pFirst.value;                   
        case Op::GreaterOrEqual: 
            return transValues.back() >= pFirst.value;                   
        default: 
            return !lp_EQUAL(transValues.back(), pFirst.value) || !lp_EQUAL(transValues[0], pFirst.value);
    }
    return false;
}

auto kernelOne = [](uint8_t t) { return t == 1; };
auto kernelZero = [](uint64_t t) { return t == 0; };
auto kernelNotZero = [](uint64_t t) { return t != 0; };

bool isTupleRangeConflict(TupleType *tupFrom, TupleType *tupTo, 
        PredIter cbegin, PredIter cend, ColumnStruct *relColumns, unsigned int pos) {
    // copy the eligible tuples into our mask vector
    vector<Metadata_t> resTuples; 
    resTuples.reserve(tupTo-tupFrom);
    resTuples.insert(resTuples.begin(), tupFrom, tupTo); 
    //resTuples.resize(tupTo-tupFrom);
    //memcpy(resTuples.data(), &*tupFrom, (tupTo-tupFrom)*sizeof(Metadata_t)); 

    size_t tplsz = relColumns[0].transactions[pos].values.size();
    vector<uint8_t> tplBitVectorRes(tplsz); uint8_t *bitvres = tplBitVectorRes.data();

    size_t activeSize = resTuples.size();
    for (; cbegin<cend; ++cbegin) {
        //cerr << "active: " << activeSize << endl;
        auto resb = resTuples.data(), rese = resTuples.data() + activeSize;

        auto& c = *cbegin;    
        auto& cTransactions = relColumns[c.column].transactions[pos];
        auto& transValues = cTransactions.values;
        //decltype(transValues.begin()) tBegin = transValues.begin(), tEnd=transValues.end();
        decltype(transValues.data()) tBegin = transValues.data(), tEnd=transValues.data()+transValues.size();
        size_t tupFromIdx{0}, tupToIdx{transValues.size()};
        switch (c.op) {
            case Op::Equal: 
                {
                    //if (transValues[0] > c.value || transValues.back() < c.value) return false;
                    auto tp = std::equal_range(tBegin, tEnd, c.value);
                    if (tp.second == tp.first) return false;
                    tupFromIdx = (tp.first - tBegin); tupToIdx = tupFromIdx + (tp.second-tp.first);
                    break;}
            case Op::Less: 
                //if (transValues[0] >= c.value) return false;
                tupToIdx -= (tEnd-std::lower_bound(tBegin, tEnd, c.value));
                if (tupToIdx == tupFromIdx) return false;
                break;
            case Op::LessOrEqual: 
                //if (transValues[0] > c.value) return false;
                tupToIdx -= (tEnd-std::upper_bound(tBegin, tEnd, c.value)); 
                if (tupToIdx == tupFromIdx) return false;
                break;
            case Op::Greater: 
                //if (transValues.back() <= c.value) return false;
                tupFromIdx += (std::upper_bound(tBegin, tEnd, c.value)-tBegin);
                if (tupToIdx == tupFromIdx) return false;
                break;
            case Op::GreaterOrEqual: 
                //if (transValues.back() < c.value) return false;
                tupFromIdx += (std::lower_bound(tBegin, tEnd, c.value)-tBegin);
                if (tupToIdx == tupFromIdx) return false;
                break;
            default: 
                // check if the active tuples have a value != to the predicate
                for (auto tpl=resb; tpl<rese; ++tpl) {
                    tpl->tuple = (tpl->tuple[c.column] == c.value) ? 0 : tpl->tuple; 
                    //tpl->tpl_id = (tpl->tuple[c.column] == c.value) ? 0 : tpl->tpl_id; 
                    //if (tpl->tuple[c.column] == c.value) tpl->tuple = 0; 
                    //if (tpl->tuple[c.column] == c.value) tpl->tpl_id = 0; 
                }
                goto LBL_CHECK_END;
        }

        // this check is done for all the operators apart from !=
        // we have to check if the active tuples are inside the result set returned
        {
            // reset the bitvector for the results and update with new results
            lp::simd::zero(bitvres, tplsz);
            auto transTuples = cTransactions.tuples.data();
            for (size_t i=tupFromIdx; i<tupToIdx; ++i) bitvres[transTuples[i].tpl_id] = (uint8_t)1;
            // remove those that are invalid
            for (auto tpl=resb; tpl<rese; ++tpl) {
                tpl->tuple = (bitvres[tpl->tpl_id]) ? tpl->tuple : 0; 
                //tpl->tpl_id = (bitvres[tpl->tpl_id]) ? tpl->tpl_id : 0; 
            }

        }
LBL_CHECK_END:
        // check if we have any valid tuple left otherwise return false
        //cerr << "active " << activeSize;
        activeSize = std::partition(resb, rese, 
                [](const Metadata_t& meta) { return meta.tuple; }) - resb;
        //[](const Metadata_t& meta) { return meta.tpl_id; }) - resb;
        //cerr << " active after " << activeSize << endl;
        //if (activeSize == 0) return false;
        //if (activeSize & (cbegin+1 == cend)) return true;
        if (activeSize < 64) return isTupleRangeConflict(resb, resb+activeSize, ++cbegin, cend);
    }
    return true;
}

static bool isTransactionConflict(const ColumnTransaction_t& transaction, Column pFirst, PredIter cbegin, PredIter cend, ColumnStruct *relColumns, unsigned int pos) {
    auto& transValues = transaction.values;
    auto& transTuples = transaction.tuples;
    decltype(transValues.data()) tBegin = transValues.data(), tEnd=transValues.data()+transValues.size();
    TupleType *tupFrom{const_cast<Metadata_t*>(transTuples.data())}, 
              *tupTo{const_cast<Metadata_t*>(transTuples.data()+transTuples.size())};
    // find the valid tuples using range binary searches based on the first predicate
    switch (pFirst.op) {
        case Op::Equal: 
            {
                //if (transValues[0] > pFirst.value || transValues.back() < pFirst.value) return false;
                auto tp = std::equal_range(tBegin, tEnd, pFirst.value);
                if (tp.second == tp.first) return false;
                tupFrom += (tp.first - tBegin); tupTo -= (tEnd-tp.second);
                break;}
        case Op::Less: 
            //if (transValues[0] >= pFirst.value) return false;
            tupTo -= (tEnd-std::lower_bound(tBegin, tEnd, pFirst.value));
            if (tupTo == tupFrom) return false;
            break;
        case Op::LessOrEqual: 
            //if (transValues[0] > pFirst.value) return false;
            tupTo -= (tEnd-std::upper_bound(tBegin, tEnd, pFirst.value)); 
            if (tupTo == tupFrom) return false;
            break;
        case Op::Greater: 
            //if (transValues.back() <= pFirst.value) return false;
            tupFrom += (std::upper_bound(tBegin, tEnd, pFirst.value)-tBegin);
            if (tupTo == tupFrom) return false;
            break;
        case Op::GreaterOrEqual: 
            //if (transValues.back() < pFirst.value) return false;
            tupFrom += (std::lower_bound(tBegin, tEnd, pFirst.value)-tBegin);
            if (tupTo == tupFrom) return false;
            break;
        default: 
            cbegin = std::prev(cbegin);
    }
    //cerr << "\n||| Transaction tuples |||" << endl;
    //cerr << "-- orig: " << transValues.size() << endl;
    //cerr << ":: after 1st: " << (tupTo-tupFrom) << endl;
    //for (auto t : transValues) cerr << t << " ";

    //cerr << "tup diff " << (tupTo - tupFrom) << endl; 
    //if (std::distance(tupFrom, tupTo) == 0) return false;

    if (tupTo - tupFrom < 64) return isTupleRangeConflict(tupFrom, tupTo, cbegin, cend);
    else return isTupleRangeConflict(tupFrom, tupTo, cbegin, cend, relColumns, pos);
    //return isTupleRangeConflict(tupFrom, tupTo, cbegin, cend);
}
*/
static bool processQueryOther(LPValidation& v, Query *q, Column *cbegin, Column *cend) {
    auto& cindex = gRelColumns[q->relationId].columns[cbegin->column].transactions;
    auto trbuckets = cindex.buckets(v.from, v.to); 
    //cerr << ri << "=" << (buckets.second - buckets.first) << endl;
    //const uint64_t rangediff = v.to - v.from;
    // find the valid tuples using range binary searches based on the first predicate
    auto csecond = cbegin + 1;

    switch (cbegin->op) {
        case Op::Less: 
            {
                for (auto cb=trbuckets.first, ce=trbuckets.second; cb!=ce; ++cb) {
                    auto tplpair = cb->lower_bound_left(cbegin->value);
                    if (tplpair.first == tplpair.second) continue;
                    if (isTupleRangeConflict(tplpair.first, tplpair.second, csecond, cend, v.from, v.to)) { return true; }
                } // end for each bucket
            }
            break;
        case Op::LessOrEqual: 
            {
                for (auto cb=trbuckets.first, ce=trbuckets.second; cb!=ce; ++cb) {
                    auto tplpair = cb->upper_bound_left(cbegin->value);
                    if (tplpair.first == tplpair.second) continue;
                    if (isTupleRangeConflict(tplpair.first, tplpair.second, csecond, cend, v.from, v.to)) { return true; }
                } // end for each bucket
            }
            break;
        case Op::Greater: 
            {
                for (auto cb=trbuckets.first, ce=trbuckets.second; cb!=ce; ++cb) {
                    auto tplpair = cb->upper_bound(cbegin->value);
                    if (tplpair.first == tplpair.second) continue;
                    if (isTupleRangeConflict(tplpair.first, tplpair.second, csecond, cend, v.from, v.to)) { return true; }
                } // end for each bucket
            }
            break;
        case Op::GreaterOrEqual: 
            {
                for (auto cb=trbuckets.first, ce=trbuckets.second; cb!=ce; ++cb) {
                    auto tplpair = cb->lower_bound(cbegin->value);
                    if (tplpair.first == tplpair.second) continue;
                    if (isTupleRangeConflict(tplpair.first, tplpair.second, csecond, cend, v.from, v.to)) { return true; }
                } // end for each bucket
            }
            break;
        case Op::NotEqual: 
            {
                for (auto cb=trbuckets.first, ce=trbuckets.second; cb!=ce; ++cb) {
                    auto tplpair = cb->lower_bound_left(cbegin->value);
                    if (isTupleRangeConflict(tplpair.first, tplpair.second, csecond, cend, v.from, v.to)) { return true; }
                    tplpair = cb->upper_bound(cbegin->value);
                    if (isTupleRangeConflict(tplpair.first, tplpair.second, csecond, cend, v.from, v.to)) { return true; }
                } // end for each bucket
            }
            break;
        case Op::Equal:
            cerr << "EQ == operator should not be executed here!!!!!" << endl;
            abort();
    }
    return false;    
}

// method that evaluates a query against any column apart from the zero-0 column 
// which is handled separately below.
static bool processQueryEQ(LPValidation& v, Query *q, Column *cbegin, Column *cend) {
    //cerr << "----- NEW VAL SESSION -----" << endl;
    auto& relColumns = gRelColumns[q->relationId].columns;
    auto& cindex = relColumns[cbegin->column].transactions;
    auto trbuckets = cindex.buckets(v.from, v.to); 
    //cerr << ri << "=" << (buckets.second - buckets.first) << endl;
    //const uint64_t rangediff = v.to - v.from;

    auto csecond = cbegin + 1;
    // TODO - cases for 1 predicate - for 2 equalities to take the less tuples
   
    // for each bucket of the column index
    for (uint32_t cbi=0, cbsz=(trbuckets.second-trbuckets.first); cbi<cbsz; ++cbi) {
        auto cb = &trbuckets.first[cbi];
        auto tplpair = cb->equal_range(cbegin->value, v.from, v.to);
        const uint64_t resdiff = tplpair.second - tplpair.first;
        // if no tuples with this value
        if (0 == resdiff) continue;
        // if result set small enough to process
        if (128 > resdiff || !(q->columnCount > 1 && csecond->op == 0)) {
            if (isTupleRangeConflict(tplpair.first, tplpair.second, csecond, cend)) { return true; }
        } else {
            // the first column bucket returned a big result
            // therefore we will try to check the second column if == predicate
            // and process the least results
            auto& cindex2 = relColumns[csecond->column].transactions;
            auto trbuckets2 = cindex2.buckets(v.from, v.to); 
            auto cb2 = &trbuckets2.first[cbi];
            
            auto tplpair2 = cb2->equal_range(csecond->value, v.from, v.to);
            const uint64_t resdiff2 = tplpair2.second - tplpair2.first;
            if (resdiff2 <= resdiff) {
                std::swap(*cbegin, *csecond);
                if (isTupleRangeConflict(tplpair2.first, tplpair2.second, csecond, cend)) { return true; }
                std::swap(*cbegin, *csecond);
            } else {
                if (isTupleRangeConflict(tplpair.first, tplpair.second, csecond, cend)) { return true; }
            }
        } // end if big result
    }    
    return false;
}

// the actual processing of a query with a predicate against primary key column
// TODO - a special index for the primary column might did a better job here!!!
static bool processQueryEQZero(LPValidation& v, Query *q, Column* cbegin, Column *cend) {
    auto& primIndex = gRelColumns[q->relationId].columns[0].transactions;
    auto trbuckets = primIndex.buckets(v.from, v.to); 
    auto csecond = cbegin + 1;
    for (auto cb=trbuckets.first, ce=trbuckets.second; cb!=ce; ++cb) {
        auto tplpair = cb->equal_range(cbegin->value);
        if (tplpair.first == tplpair.second) continue;
        if (isTupleRangeConflict(tplpair.first, tplpair.second, csecond, cend, v.from, v.to)) { return true; }
    } // for all buckets
    return false;
}

// this method is responsible to process all the equality queries of column CI
// in relation RI.
static void processEqualityQueries(uint32_t tid, uint32_t ri, uint32_t ci) { (void)tid;
    auto& rq = gRelQ[ri].columns[ci].queries;
    if (rq.empty()) return;
    //cerr << "rel: " << ri << " col " << ci << " ==: " << rq.size() << endl;
    QMeta_t *qb = rq.data(), *qe = rq.data()+rq.size();

    auto processFunc = (ci!=0) ? &processQueryEQ : &processQueryEQZero;

    for (; qb!=qe; ) {
        auto& cmeta = *qb++;
        const size_t resPos = cmeta.lpv->validationId - resIndexOffset;
        if (gPendingResults[resPos]) { continue; }
        if (processFunc(*cmeta.lpv, cmeta.rq, (Column*)cmeta.rq->columns, 
                    ((Column*)cmeta.rq->columns)+cmeta.rq->columnCount)){
            gPendingResults[resPos] = true;
        }
    } // end of all queries
}


static bool isValidationConflict(LPValidation& v) {
    //cerr << "\n========= validation " << v.validationId << " =========" << endl; 
    //cerr << "qc: " << vq.queryCount << " from: " << vq.from << " to: " << vq.to << endl; 
    for (auto q : v.queries) {
        Query& rq = *q;
        if (unlikely(rq.columnCount == 0)) { 
            //cerr << "empty: " << v.validationId << endl; 
            auto& transactionsCheck = gRelations[rq.relationId].transLogTuples;
            auto transFromCheck = std::lower_bound(transactionsCheck.begin(), 
                    transactionsCheck.end(), v.from, TransLogComp);
            if (transFromCheck == transactionsCheck.end() || transFromCheck->first > v.to) {
                // no transactions exist for this query
                continue;
            } else { 
                // transactions exist for this query so it is a conflict
                return true;
            }; 
        }
        uint32_t colCountUniq = rq.columnCount; 
        auto cbegin = reinterpret_cast<Query::Column*>(rq.columns),
             cend = cbegin + colCountUniq;
        if (processQueryOther(v, &rq, cbegin, cend)) { return true; }
        else continue;
/*
        // --- PROCESSING BEFORE THE BUCKETTED INDEX --- 2 DAYS before the end

        auto pFirst = *reinterpret_cast<Query::Column*>(rq.columns);
        // just find the range of transactions we want in this relation
        auto& relColumns = gRelColumns[rq.relationId].columns;
        auto& transactions = relColumns[pFirst.column].transactions;
        auto transFrom = std::lower_bound(transactions.begin(), transactions.end(), v.from, CTRSLessThan);
        auto transTo = std::upper_bound(transFrom, transactions.end(), v.to, CTRSLessThan);

        //cerr << (transTo - transFrom) << endl;

        uint32_t pos = std::distance(transactions.begin(), transFrom);

        // increase cbegin to point to the 2nd predicate to avoid the increment inside the function
        auto cbSecond = cbegin+1;

           if (colCountUniq > 2) {
           auto& cb=cbegin[0], cb1=cbegin[1], cb2=cbegin[2];
           for(; transFrom<transTo; ++transFrom, ++pos) {  
           if (    !((!(cb.op) & !lp_EQUAL((relColumns[cb.column].transactionsORs[pos] & cb.value), cb.value))
           || (!(cb1.op) & !lp_EQUAL((relColumns[cb1.column].transactionsORs[pos] & cb1.value), cb1.value))
           || (!(cb2.op) & !lp_EQUAL((relColumns[cb2.column].transactionsORs[pos] & cb2.value), cb2.value)))
           && isTransactionConflict(*transFrom, pFirst, cbSecond, cend, relColumns.get(), pos)) { return true; }
           } // end of all the transactions for this relation for this specific query
           } else if (colCountUniq > 1) {
               auto& cb=cbegin[0], cb1=cbegin[1];
               for(; transFrom<transTo; ++transFrom, ++pos) {  
                   if (    !((!(cb.op) && !lp_EQUAL((relColumns[cb.column].transactionsORs[pos] & cb.value), cb.value))
                               || (!(cb1.op) && !lp_EQUAL((relColumns[cb1.column].transactionsORs[pos] & cb1.value), cb1.value)))
                           && isTransactionConflict(*transFrom, pFirst, cbSecond, cend, relColumns.get(), pos)) { return true; }
               } // end of all the transactions for this relation for this specific query
           } else  {
               auto& cb = cbegin[0];
               if (!cb.op) { 
                   for(; transFrom<transTo; ++transFrom, ++pos) {  
                       if (!(!lp_EQUAL((relColumns[cb.column].transactionsORs[pos] & cb.value), cb.value))
                               && isTransactionConflict(*transFrom, pFirst)) { return true; }
                   } // end of all the transactions for this relation for this specific query
               } else {
                   for(; transFrom<transTo; ++transFrom) {  
                       if (isTransactionConflict(*transFrom, pFirst)) { return true; }
                   } // end of all the transactions for this relation for this specific query
               }
           }
*/
        //for(; transFrom<transTo; ++transFrom, ++pos) {  
        //    if (isTransactionConflict(*transFrom, pFirst, cbSecond, cend, relColumns.get(), pos)) { return true; }
        //} // end of all the transactions for this relation for this specific query
     
    }// end for all queries
    return false;
}

void processEqualityQ(uint32_t nThreads, uint32_t tid, void *args) {
    (void)tid; (void)nThreads; (void)args;// to avoid unused warning
#ifdef LPDEBUG
    auto qproc = LPTimer.getChrono();
#endif
    uint64_t totalCols = gEQCols.size();
    uint32_t ri, col;
    for (uint64_t rc = gNextEQCol++; rc < totalCols; rc=gNextEQCol++) {
        lp::validation::unpackRelCol(gEQCols[rc], ri, col);
        processEqualityQueries(tid, ri, col);
    } // end of while columns to update     
#ifdef LPDEBUG
    LPTimer.validationsProcessingIndex += LPTimer.getChrono(qproc);
#endif
}

/**
  * MAIN task executed by threadpool
  */
void processPendingValidationsTask(uint32_t nThreads, uint32_t tid, void *args) {
    (void)tid; (void)nThreads; (void)args;// to avoid unused warning

    // process all the queries with == predicates
    processEqualityQ(nThreads, tid, args);

    uint64_t totalPending = gPendingValidations.size();
    // get a validation ID - atomic operation
    for (uint64_t vi = gNextPending++; vi < totalPending; vi=gNextPending++) {
        auto& v = gPendingValidations[vi];
        uint64_t resPos = v.validationId - resIndexOffset;
        auto& atoRes = gPendingResults[resPos];
        if (atoRes) continue;
        if(isValidationConflict(v)) { atoRes = true; }
    } // while true take more validations 
}

