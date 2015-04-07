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
//#include "include/MultiTaskPool.hpp"

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
/*
   std::ostream& operator<< (std::ostream& os, const LPQuery& o) {
   os << "{" << o.relationId << "-" << o.columnCount << ":: " << o.predicates << "::" << o.satisfiable << "}";
   return os;
   }
 */
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

template<typename T>
struct LPPair {
    T *a;
    T *b;

    LPPair() : a(nullptr), b(nullptr) {}
    LPPair(T *v) : a(v), b(nullptr) {}
    LPPair(T *v, T *v2) : a(v), b(v2) {}

    void put(T *v) {
        if (a == nullptr) a = v;
        else b = v;
    }
};
template<>
struct LPPair<uint64_t> {
    uint64_t a;
    tuple_t tpl_a;
    uint64_t b;
    tuple_t tpl_b;
    uint64_t c;
    tuple_t tpl_c;

    LPPair() : a(UINT64_MAX), tpl_a(nullptr), b(UINT64_MAX), tpl_b(nullptr), c(UINT64_MAX), tpl_c(nullptr) {}
    //LPPair(uint64_t v) : a(v), b(UINT64_MAX), tpl_a(nullptr), tpl_b(nullptr) {}
    //LPPair(uint64_t v, uint64_t v2) : a(v), b(v2), tpl_a(nullptr), tpl_b(nullptr) {}

    void put_u(uint64_t v, tuple_t _tpl) {
        if (a == UINT64_MAX) { a=v; tpl_a=_tpl; }
        else if (b == UINT64_MAX) { b=v; tpl_b=_tpl; }
        else { c=v; tpl_c=_tpl; }
    }
};

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
    CIndex transactions;
    //vector<uint64_t> transactionsORs;
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
    ColumnStruct *columns;
};
static RelationColumns *gRelColumns;
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
    CQ_t *columns;
};
static RelQ_t *gRelQ;

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
static RelationStruct *gRelations;

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
static RelTransLock *gRelTransMutex;
static vector<TRMapPhase> *gTransParseMapPhase;

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
//static inline void checkPendingTransactions(ISingleTaskPool *pool);
void processPendingIndexTask(uint32_t nThreads, uint32_t tid, void *args);

///--------------------------------------------------------------------------
///--------------------------------------------------------------------------

static void processDefineSchema(const DefineSchema& d) {
    gSchema.reset(new uint32_t[d.relationCount]);
    memcpy(gSchema.get(), d.columnCounts, sizeof(uint32_t)*d.relationCount);
    NUM_RELATIONS = d.relationCount;
    //cerr << "relations: " << NUM_RELATIONS << endl;

    gRelTransMutex = new RelTransLock[d.relationCount];
    gTransParseMapPhase = new vector<TRMapPhase>[d.relationCount];

    gRelations = new RelationStruct[d.relationCount];
    gRelColumns = new RelationColumns[d.relationCount];
    gRelQ = new RelQ_t[d.relationCount];
    //cerr << endl << "relations: " << NUM_RELATIONS << endl;
    const uint32_t rels = d.relationCount;
    for(uint32_t ri=0; ri<rels; ++ri) {
        //cerr << " " << gSchema[ci];
        //gRelColumns[ri].columns.reset(new ColumnStruct[gSchema[ri]]);
        gRelColumns[ri].columns = new ColumnStruct[gSchema[ri]];
        gRelQ[ri].columns = new CQ_t[gSchema[ri]];
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

    // try to put all the queries into a vector
    /*    
          const char* qreader=v.queries;
    //cerr << "\n----- val: " << v.validationId << " : " << v.from << "-" << v.to << " ------" << endl;
    for (uint32_t i=0;i<v.queryCount;++i) {
    const Query *q=reinterpret_cast<const Query*>(qreader);
    //cerr << "\t| q: " << q->relationId << endl;
    for (size_t ci=0; ci<q->columnCount; ++ci) {
    //cerr << q->columns[ci] << " ";
    cerr << q->columns[ci].column << " :" << q->columns[ci].op << endl;
    }
    //cerr << endl;
    qreader+=sizeof(Query)+(sizeof(Query::Column)*q->columnCount);
    }
     */

    //cerr << v.validationId << "====" << v.from << ":" << v.to << "=" << v.queryCount << endl;
    //gPendingValidationsMutex.lock();
    gPendingValidations.emplace_back(v.validationId, v.from, v.to, v.queryCount, msg);    
    // update the global pending validations to reflect this new one
    ++gPVunique;
    //gPendingValidationsMutex.unlock();

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
            //cerr << (vp.second == true ? one : zero); 
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

    // clean the index columns
    /*
    auto ub = upper_bound(cRelCol.columns[0].transactions.begin(), cRelCol.columns[0].transactions.end(),
            trans_id,
            [](const uint64_t target, const ColumnTransaction_t& ct){ return target < ct.trans_id; });
    size_t upto = std::distance(cRelCol.columns[0].transactions.begin(), ub);
    */
        
    //cerr << ri << " = " << cRelCol.columns[0].transactions.size() << " | " << cRelCol.columns[1].transactions.size() << endl;
    for (uint32_t ci=0,sz=gSchema[ri]; ci<sz; ++ci) {
        cRelCol.columns[ci].transactions.erase(trans_id);
    }
    
    /*
    // clean the transactions log
    auto& transLog = gRelations[i].transLog;         
    //cerr << "size bef: " << transLog.size() << endl;
    for (auto it = transLog.begin(), tend=transLog.end(); it!=tend && ((*it)->trans_id <= f.transactionId); ) {
    if ((*it)->aliveTuples == 0 && (*it)->last_del_id <= f.transactionId) { it = transLog.erase(it); tend=transLog.end(); }
    else ++it;
    }
     */
}

static atomic<uint64_t> gNextFRel;
void processForgetThreaded(uint32_t nThreads, uint32_t tid, void *args) { 
    (void)tid; (void)nThreads; (void)args;// to avoid unused warning
    auto f = *reinterpret_cast<Forget*>(args);
    for (uint64_t ri = gNextFRel++; ri < NUM_RELATIONS; ri=gNextFRel++) { 
        forgetRel(f.transactionId, ri);
    }
}

static void processForget(const Forget& f, ISingleTaskPool* pool) {
#ifdef LPDEBUG
    auto start = LPTimer.getChrono();
#endif
    /*
       gNextFRel = 0; 
       pool->startSingleAll(processForgetThreaded, (void*)&f);
       pool->waitSingleAll();
     */

    (void)pool;
    // delete the transactions from the columns index
    for (uint32_t i=0; i<NUM_RELATIONS; ++i) {
        forgetRel(f.transactionId, i);
    }

#ifdef LPDEBUG
    LPTimer.forgets += LPTimer.getChrono(start);
#endif
}
//---------------------------------------------------------------------------





struct FR_t{
    //std::atomic<bool> f1;
    //std::atomic<bool> f2;
    //std::atomic<uint64_t> colsdone;
    //std::atomic<uint64_t> cols;
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
std::mutex mMtxF1;
std::condition_variable mCondF1;
std::atomic<uint32_t> finishedF1;

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
        //msgReader = ReaderIOFactory::createAsync(stdin);
        msgReader = ReaderIOFactory::create(stdin);
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
    //SingleTaskPool workerThreads2(numOfThreads>>1, processPendingValidationsTask);
    //SingleTaskPool workerThreads2(1, processPendingValidationsTask);
    //workerThreads2.initThreads();

    cerr << "ColumnStruct: " << sizeof(ColumnStruct) << " RelTransLog: " << sizeof(RelTransLog) << " RelationStruct: " << sizeof(RelationStruct) << " CTransStruct: " << sizeof(CTransStruct) << endl;

    try {

#ifdef LPDEBUG
        uint64_t totalForgets = 0, totalFlushes = 0;
#endif

        while (true) {
#ifdef LPDEBUG // I put the inner timer here to avoid stalls in the msgQ
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
                        //workerThreads2.destroy();
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
        // TODO - lock here to make it to make all the deletions parallel naive locking first - 
        // TODO try to lock with try_lock and try again at the end if some relations failed
        {// start of lock_guard
            uint64_t *ptr = new uint64_t[o.rowCount];
            for (uint32_t c=0; c<o.rowCount; ++c) ptr[c] = o.keys[c];
            //memcpy(ptr, o.keys, sizeof(uint64_t)*o.rowCount);
            //gRelTransMutex[o.relationId].lock(); 
            gTransParseMapPhase[o.relationId].emplace_back(t.transactionId, true, o.rowCount, ptr);
            //gRelTransMutex[o.relationId].unlock(); 
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
        // TODO - lock here to make it to make all the deletions parallel naive locking first - 
        // TODO try to lock with try_lock and try again at the end if some relations failed
        {// start of lock_guard
            uint64_t* tptr = new uint64_t[relCols*o.rowCount];
            uint32_t sz = relCols*o.rowCount;
            for (uint32_t c=0; c<sz; ++c) tptr[c] = o.values[c];
            //memcpy(tptr, o.values, sizeof(uint64_t)*relCols*o.rowCount);
            //gRelTransMutex[o.relationId].lock(); 
            gTransParseMapPhase[o.relationId].emplace_back(t.transactionId, false, o.rowCount, tptr);
            //gRelTransMutex[o.relationId].unlock(); 
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


void processUpdateIndexTask(uint32_t nThreads, uint32_t tid, void *args);
static void updateRelCol(uint32_t tid, uint32_t ri, uint32_t col);

/*
   struct TRMapPhase {
   uint64_t trans_id;
   bool isDelOp;
   uint32_t rowCount;
   uint64_t *values; // delete op => row keys to delete | insert => tuples
   }
//static unique_ptr<vector<TRMapPhase>[]> gTransParseMapPhase;
 */
static std::atomic<uint64_t> gNextIndex;  // for the index
static std::atomic<uint32_t> gNextReqCol; // for the update columns

void processPendingIndexTask(uint32_t nThreads, uint32_t tid, void *args) {
    (void)tid; (void)nThreads; (void)args;// to avoid unused warning
    //cerr << "::: tid " << tid << "new" << endl;

    for (uint64_t ri = gNextIndex++; ri < NUM_RELATIONS; ri=gNextIndex++) {

        // take the vector with the transactions and sort it by transaction id in order to apply them in order
        auto& relTrans = gTransParseMapPhase[ri];
        if (unlikely(relTrans.empty())) { 
            //mMtxF1.lock();
            gFR[ri].f1 = true; gFR[ri].f2 = true;
            ++finishedF1; 
            //mCondF1.notify_all();
            //mMtxF1.unlock();
            continue; 
        }

        //cerr << "tid " << tid << " got " << ri << " = " << relTrans.size() << endl;

        //std::sort(relTrans.begin(), relTrans.end(), TRMapPhaseByTrans);

        auto& relation = gRelations[ri];
        uint32_t relCols = gSchema[ri];

        uint64_t lastTransId = relTrans[0].trans_id;

        //auto& primIndex = relation.primaryIndex;
        //auto firstTrans = relTrans[0].trans_id;
        //CIndex::Bucket *trb=primIndex.bucketNext(firstTrans);

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

                //trb=primIndex.bucketNext(trans.trans_id);
            }

            if (trans.isDelOp) {
                // this is a delete operation
                for (const uint64_t* key=trans.values,*keyLimit=key+trans.rowCount;key!=keyLimit;++key) {
                    auto lb = relation.insertedRows.find(*key);
                    if (lb != relation.insertedRows.end()) {
                        // lb->second is a pair<uint64_t, uint64_t*> - trans_id/tuple
                        // decrease counter of trans tuples
                        //auto tit = lower_bound(relation.transLog.begin(), relation.transLog.end(), lb->second.first, RTLComp);
                        //(*tit)->last_del_id = trans.trans_id;
                        //--(*tit)->aliveTuples;

                        // update the relation transactions - transfer ownership of the tuple
                        //tuple_t tpl = lb->second.second;
                        //operations.push_back(tpl);
                        operations.push_back(lb->second.second);

                        // remove the row from the relations table 
                        relation.insertedRows.erase(lb);

                        // insert the value into the primary key index
                        //trb->insert(trans.trans_id, operations.back(), operations.back()[0]);
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

                    // insert the value into the primary key index
                    //trb->insert(trans.trans_id, vals, vals[0]);
                }
                // TODO - THIS HAS TO BE IN ORDER - each relation will have its own transaction history from now on
                relation.transLog.emplace_back(new RelTransLog(trans.trans_id, trans.values, trans.rowCount));
            }
        }
        // store the last transaction data
        // store the operations for the last transaction
        if (likely(!operations.empty()))
            relation.transLogTuples.emplace_back(lastTransId, move(operations));

        // sort the index - HERE CAN DO - inplace merge instead
        //primIndex.sortFrom(firstTrans);

        //mMtxF1.lock();
        gFR[ri].f1 = true;
        ++finishedF1; 
        //mCondF1.notify_all();
        //mMtxF1.unlock();
    } // end of all relations

    // TODO - MERGE THE UPDATE INDEX - WITH UPDATE COLUMNS - TODO
    //processUpdateIndexTask(nThreads, tid, nullptr);
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
/*
    //if (col == 0) {
    auto trp = transFrom;
    // while we haven't processed all transactions
    while (trp<tEnd) {
        // will return the start of the bucket I should write and the last transaction I should stop
        auto bpair = cindex.prepareBucket(col==0, trp, tEnd);
        auto tplPtr = bpair.first;
        auto currentEnd = bpair.second;
        for (; trp<currentEnd; ++trp) {
            auto trans_id = trp->first;
            for (auto tpl : trp->second) { 
                *tplPtr++ = {tpl[col], trans_id, tpl}; 
            }
        }
    } // while still transactions to process
    } else {*/
    // for all the transactions in the relation
    for(auto trp=transFrom; trp<tEnd; ++trp) {
        // allocate vectors for the current new transaction to put its data
        auto trans_id = trp->first;
        CIndex::Bucket &trb = *cindex.bucketNext(trans_id, col==0);
       
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
/*
void processUpdateIndexTask(uint32_t nThreads, uint32_t tid, void *args) {
    (void)tid; (void)nThreads; (void)args;// to avoid unused warning
    uint64_t totalCols = gRequiredColumns.size();
    for (uint64_t rc = gNextReqCol++; rc < totalCols; rc=gNextReqCol++) {
        uint32_t ri, col;
        lp::validation::unpackRelCol(gRequiredColumns[rc], ri, col);
        updateRelCol(tid, ri, col);
    } // end of while columns to update     
}

static inline void checkPendingTransactions(ISingleTaskPool *pool) { (void)pool;
#ifdef LPDEBUG
    auto startIndex = LPTimer.getChrono();
#endif
    //cerr << "::: session start ::::" << endl;
    gNextIndex = 0;
    //processPendingIndexTask(1,0,nullptr); 
    pool->startSingleAll(processPendingIndexTask);
    pool->waitSingleAll();

    for (uint32_t r=0; r<NUM_RELATIONS; ++r) gTransParseMapPhase[r].clear();
#ifdef LPDEBUG
    LPTimer.transactionsIndex += LPTimer.getChrono(startIndex);
#endif

#ifdef LPDEBUG
    auto startUpdIndex = LPTimer.getChrono();
#endif
    gNextReqCol = 0;
    //processUpdateIndexTask(1, 0, nullptr);
    pool->startSingleAll(processUpdateIndexTask);
    pool->waitSingleAll();
#ifdef LPDEBUG
    LPTimer.updateIndex += LPTimer.getChrono(startUpdIndex);
#endif
}
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



void processUpdateIndexTask(uint32_t nThreads, uint32_t tid, void *args) {
    (void)tid; (void)nThreads; (void)args;// to avoid unused warning
    //uint64_t totalCols = gRequiredColumns.size();
    //for (uint64_t rc = gNextReqCol++; rc < totalCols; rc=gNextReqCol++) {
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
        
        //mMtxF1.lock();
        //if (finishedF1 < NUM_RELATIONS) mCondF1.wait(mMtxF1, []{return true;});
        if (finishedF1 < NUM_RELATIONS) lp_spin_sleep(std::chrono::microseconds(10));
        else { allfinished = true; }
        //mMtxF1.unlock();
    } // all relations finished
}

void parallelTask1(uint32_t nThreads, uint32_t tid, void *args) { (void)nThreads; (void)tid; (void)args;
    processPendingIndexTask(nThreads, tid, args); // F1
    if (tid < 4) {
        createQueryIndexTask(nThreads, tid, args);    // F3
    }
    processUpdateIndexTask(nThreads, tid, args);
}


static void checkPendingValidations(ISingleTaskPool *pool, ISingleTaskPool *pool2) {
    if (unlikely(gPendingValidations.empty())) return;
    if (unlikely(gFR.empty())) gFR.resize(NUM_RELATIONS);
#ifdef LPDEBUG
    auto startQuery = LPTimer.getChrono();
#endif
    /*
    gNextPending = 0;
    pool2->startSingleAll(createQueryIndexTask);
    // check if there is any pending index creation to be made before checking validation
    checkPendingTransactions(pool);
    finishQueryIndex(pool2);
    */

    
    resetUpdateStats();
    // trans-index & qindex
    gNextPending = 0; // F3
    gNextIndex = 0; // F1
    gNextReqCol = 0; // - F2 
    pool->startSingleAll(parallelTask1); // F1 + F3
    pool->waitSingleAll(); 
    //finishQueryIndex(pool); // F3

    //pool->startSingleAll(processUpdateIndexTask); // F2

    // this is needed by the trans-index after it has finished - F1
    for (uint32_t r=0; r<NUM_RELATIONS; ++r) gTransParseMapPhase[r].clear();
    
    //pool->waitSingleAll(); // F2
    













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
    //memset(gPendingResults.data(), 0, sizeof(PendingResultType)*gPRsz);
    for (auto gpr=gPendingResults.data(), gpre=gpr+gPRsz; gpr<gpre; ) *gpr++ = 0;

    gNextPending = 0; gNextEQCol = 0;
    //processPendingValidationsTask(1,0,nullptr);
    pool->startSingleAll(processPendingValidationsTask);
    //pool2->startSingleAll(processPendingValidationsTask);
    //pool2->waitSingleAll();
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

//typedef tuple_t TupleType;
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

bool ALWAYS_INLINE isTupleRangeConflict(TupleType *tupFrom, TupleType *tupTo, PredIter cbegin, PredIter cend) {
    if (cbegin == cend && tupTo-tupFrom != 0) return true;
    for(; tupFrom!=tupTo; ++tupFrom) {
        if (isTupleConflict(cbegin, cend, tupFrom->tuple)) return true;
    } // end of all tuples for this transaction
    return false;
    //return std::find_if(tupFrom, tupTo, [=](TupleType& tup) { return isTupleConflict(cbegin, cend, tup);}) != tupTo;
}

/*
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
    const uint64_t rangediff = v.to - v.from;
    // find the valid tuples using range binary searches based on the first predicate
    auto csecond = cbegin + 1;

    switch (cbegin->op) {
        case Op::Less: 
            {
                for (auto cb=trbuckets.first, ce=trbuckets.second; cb!=ce; ++cb) {
                    auto tplpair = cb->lower_bound_left(cbegin->value);
                    if (tplpair.first == tplpair.second) continue;
                    for (auto ctpl=tplpair.first, tple=tplpair.second; ctpl<tple; ++ctpl) {
                        if ( (uint64_t)(ctpl->trans_id - v.from) <= (rangediff)) {
                            if (isTupleConflict(csecond, cend, ctpl->tuple)) {
                                return true;;
                            }
                        }
                    } // end for all trans in the bucket
                } // end for each bucket
            }
            //tupTo -= (tEnd-std::lower_bound(tBegin, tEnd, pFirst.value));
            //if (tupTo == tupFrom) return false;
            break;
        case Op::LessOrEqual: 
            {
                for (auto cb=trbuckets.first, ce=trbuckets.second; cb!=ce; ++cb) {
                    auto tplpair = cb->upper_bound_left(cbegin->value);
                    if (tplpair.first == tplpair.second) continue;
                    for (auto ctpl=tplpair.first, tple=tplpair.second; ctpl<tple; ++ctpl) {
                        if ( (uint64_t)(ctpl->trans_id - v.from) <= (rangediff)) {
                            if (isTupleConflict(csecond, cend, ctpl->tuple)) {
                                return true;;
                            }
                        }
                    } // end for all trans in the bucket
                } // end for each bucket
            }
            //tupTo -= (tEnd-std::upper_bound(tBegin, tEnd, pFirst.value)); 
            //if (tupTo == tupFrom) return false;
            break;
        case Op::Greater: 
            {
                for (auto cb=trbuckets.first, ce=trbuckets.second; cb!=ce; ++cb) {
                    auto tplpair = cb->upper_bound(cbegin->value);
                    if (tplpair.first == tplpair.second) continue;
                    for (auto ctpl=tplpair.first, tple=tplpair.second; ctpl<tple; ++ctpl) {
                        if ( (uint64_t)(ctpl->trans_id - v.from) <= (rangediff)) {
                            if (isTupleConflict(csecond, cend, ctpl->tuple)) {
                                return true;;
                            }
                        }
                    } // end for all trans in the bucket
                } // end for each bucket
            }
            //tupFrom += (std::upper_bound(tBegin, tEnd, pFirst.value)-tBegin);
            //if (tupTo == tupFrom) return false;
            break;
        case Op::GreaterOrEqual: 
            {
                for (auto cb=trbuckets.first, ce=trbuckets.second; cb!=ce; ++cb) {
                    auto tplpair = cb->lower_bound(cbegin->value);
                    if (tplpair.first == tplpair.second) continue;
                    for (auto ctpl=tplpair.first, tple=tplpair.second; ctpl<tple; ++ctpl) {
                        if ( (uint64_t)(ctpl->trans_id - v.from) <= (rangediff)) {
                            if (isTupleConflict(csecond, cend, ctpl->tuple)) {
                                return true;;
                            }
                        }
                    } // end for all trans in the bucket
                } // end for each bucket
            }
            //tupFrom += (std::lower_bound(tBegin, tEnd, pFirst.value)-tBegin);
            //if (tupTo == tupFrom) return false;
            break;
        case Op::NotEqual: 
            {
                for (auto cb=trbuckets.first, ce=trbuckets.second; cb!=ce; ++cb) {
                    auto tplpair = cb->lower_bound_left(cbegin->value);
                    for (auto ctpl=tplpair.first, tple=tplpair.second; ctpl<tple; ++ctpl) {
                        if ( (uint64_t)(ctpl->trans_id - v.from) <= (rangediff)) {
                            if (isTupleConflict(csecond, cend, ctpl->tuple)) {
                                return true;;
                            }
                        }
                    } // end for all trans in the bucket
                    tplpair = cb->upper_bound(cbegin->value);
                    for (auto ctpl=tplpair.first, tple=tplpair.second; ctpl<tple; ++ctpl) {
                        if ( (uint64_t)(ctpl->trans_id - v.from) <= (rangediff)) {
                            if (isTupleConflict(csecond, cend, ctpl->tuple)) {
                                return true;;
                            }
                        }
                    } // end for all trans in the bucket
                } // end for each bucket
            }
            break;
        case Op::Equal:
            cerr << "EQ == operator should not be executed here!!!!!" << endl;
            abort();
    }
    //cerr << "\n||| Transaction tuples |||" << endl;
    //cerr << "-- orig: " << transValues.size() << endl;
    //cerr << ":: after 1st: " << (tupTo-tupFrom) << endl;
    //for (auto t : transValues) cerr << t << " ";
    return false;    
}

static bool processQueryEQ(LPValidation& v, Query *q, Column *cbegin, Column *cend) {
    //cerr << "----- NEW VAL SESSION -----" << endl;
    //uint64_t cnt = 0, cnt2 = 0;
    auto& relColumns = gRelColumns[q->relationId].columns;
    auto& cindex = relColumns[cbegin->column].transactions;
    auto trbuckets = cindex.buckets(v.from, v.to); 
    //cerr << ri << "=" << (buckets.second - buckets.first) << endl;
    //const uint64_t rangediff = v.to - v.from;

    auto csecond = cbegin + 1;
    // TODO - cases for 1 predicate - for 2 equalities to take the less tuples
    
    if (q->relationId == 3 && cbegin->column == 4 && q->columnCount > 1 && csecond->op == 0) {
        std::swap(*cbegin, *csecond);
        auto& cindex = relColumns[cbegin->column].transactions;
        auto trbuckets = cindex.buckets(v.from, v.to); 
    for (auto cb=trbuckets.first, ce=trbuckets.second; cb!=ce; ++cb) {
        auto tplpair = cb->equal_range(cbegin->value, v.from, v.to);
        if (tplpair.first == tplpair.second) continue;
        //cerr<< "found : " << (rp.second-rp.first) << endl;
        //register const size_t tplsz=tplpair.second-tplpair.first;
        auto ctpl = tplpair.first;
        // know that all the tuples we got are in the range we want - equal_range in bucket guarantees that
        for (const auto tple=tplpair.second; (ctpl < tple); ) { 
            if (isTupleConflict(csecond, cend, (ctpl++)->tuple)) { return true; }
        }
        //cerr<< "break: " << (tplpair.second-ctpl) << "/" << (tplpair.second-tplpair.first) << endl;  
    } // for all buckets
    } else {
    for (auto cb=trbuckets.first, ce=trbuckets.second; cb!=ce; ++cb) {
        auto tplpair = cb->equal_range(cbegin->value, v.from, v.to);
        if (tplpair.first == tplpair.second) continue;
        //cerr<< "found : " << (rp.second-rp.first) << endl;
        //register const size_t tplsz=tplpair.second-tplpair.first;
        auto ctpl = tplpair.first;
        // know that all the tuples we got are in the range we want - equal_range in bucket guarantees that
        for (const auto tple=tplpair.second; (ctpl < tple); ) { 
            if (isTupleConflict(csecond, cend, (ctpl++)->tuple)) { return true; }
        }
        //cerr<< "break: " << (tplpair.second-ctpl) << "/" << (tplpair.second-tplpair.first) << endl;  
    } // for all buckets
    }
    return false;
}
static bool processQueryEQZero(LPValidation& v, Query *q, Column* cbegin, Column *cend) {
    //cerr << "----- NEW VAL SESSION -----" << endl;
    auto& primIndex = gRelColumns[q->relationId].columns[0].transactions;
    //cerr << ri << "=" << (buckets.second - buckets.first) << endl;
    auto trbuckets = primIndex.buckets(v.from, v.to); 
    const uint64_t rangediff = v.to - v.from;
    
    auto csecond = cbegin + 1;

    //cerr << "query: " << cmeta.from << "-" << cmeta.to <<endl;
    //cerr << "brange " << trp.first->trmin << "-" << trp.first->trmax << " & " << trp.second->trmin << "-" << trp.second->trmax << endl;
    for (auto cb=trbuckets.first, ce=trbuckets.second; cb!=ce; ++cb) {
        auto tplpair = cb->equal_range(cbegin->value);
        if (tplpair.first == tplpair.second) continue;
        for (auto ctpl=tplpair.first, tple=tplpair.second; ctpl<tple; ++ctpl) {
            //if (ctpl->trans_id <= cmeta.to && ctpl->trans_id >= cmeta.from) {
            if ( (uint64_t)(ctpl->trans_id - v.from) <= (rangediff)) {
                if (isTupleConflict(csecond, cend, ctpl->tuple)) {
                    return true;
                }
            }
        }
    } // for all buckets

    return false;
}

static void processEqualityQueries(uint32_t tid, uint32_t ri, uint32_t ci) {
    auto& rq = gRelQ[ri].columns[ci].queries;
    if (rq.empty()) return;
    //cerr << "rel: " << ri << " col " << ci << " ==: " << rq.size() << endl;
    QMeta_t *qb = rq.data(), *qe = rq.data()+rq.size();

    auto processFunc = (ci!=0) ? &processQueryEQ : &processQueryEQZero;

    for (; qb!=qe; ) {
        auto& cmeta = *qb++;
        const size_t resPos = cmeta.lpv->validationId - resIndexOffset;
        if (gPendingResults[resPos]) { continue; }
        if (processFunc(*cmeta.lpv, cmeta.rq, (Column*)cmeta.rq->columns, ((Column*)cmeta.rq->columns)+cmeta.rq->columnCount)){
            gPendingResults[resPos] = true;
        }
    } // end of all queries
}



static bool isValidationConflict(LPValidation& v) {
    // TODO - MAKE A PROCESSING OF THE QUERIES AND PRUNE SOME OF THEM OUT
    /*
       cerr << "\n========= validation " << v.validationId << " =========" << endl; 
       cerr << "qc: " << vq.queryCount << " from: " << vq.from << " to: " << vq.to << endl; 
     */

    for (auto q : v.queries) {
        Query& rq = *q;
    //const ValidationQueries& vq = *reinterpret_cast<ValidationQueries*>(v.rawMsg->data.data());
    //const char* qreader = vq.queries; uint32_t columnCount;
    //for (uint32_t i=0; i<vq.queryCount; ++i, qreader+=sizeof(Query)+(sizeof(Query::Column)*columnCount)) {
        //Query& rq=*const_cast<Query*>(reinterpret_cast<const Query*>(qreader));
        //columnCount = rq.columnCount;
        //cerr << " " << i;

        if (unlikely(rq.columnCount == 0)) { 
            //cerr << "empty: " << v.validationId << endl; 
            auto& transactionsCheck = gRelations[rq.relationId].transLogTuples;
            auto transFromCheck = std::lower_bound(transactionsCheck.begin(), transactionsCheck.end(), v.from, TransLogComp);
            if (transFromCheck == transactionsCheck.end() || transFromCheck->first > v.to) {
                // no transactions exist for this query
                continue;
            } else { 
                // transactions exist for this query so it is a conflict
                return true;
            }; 
        }

#ifdef LPDEBUG
//auto startInner = LPTimer.getChrono();
#endif 
        //uint32_t colCountUniq = lp::query::preprocess(rq); 
        //if (!lp::query::satisfiable(&rq, colCountUniq)) { continue; } // go to the next query
        //if (!lp::query::preprocess(rq, gSchema[rq.relationId])) { continue; }
#ifdef LPDEBUG
//LPTimer.satCheck += LPTimer.getChrono(startInner);
#endif 
    
        uint32_t colCountUniq = rq.columnCount; 
    
        auto cbegin = reinterpret_cast<Query::Column*>(rq.columns),
             cend = cbegin + colCountUniq;
        
        /*if (cbegin->op == Op::Equal) {
            if (cbegin->column == 0) {
                if (processQueryEQZero(v, &rq, cbegin, cend)) { return true; }  
                else continue;
            } else {
                if (processQueryEQ(v, &rq, cbegin, cend)) { return true; }  
                else continue;
            }
        } else {*/
            if (processQueryOther(v, &rq, cbegin, cend)) { return true; }
            else continue;
        //}
    
        /*
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
        */

        /*
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
        //cerr << ":: val " << v.validationId << endl;
    
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
    //uint64_t totalCols = gRequiredColumns.size();
    uint64_t totalCols = gEQCols.size();
    uint32_t ri, col;
    for (uint64_t rc = gNextEQCol++; rc < totalCols; rc=gNextEQCol++) {
        lp::validation::unpackRelCol(gEQCols[rc], ri, col);
        processEqualityQueries(tid, ri, col);
    } // end of while columns to update     
    /*
    for (uint32_t ri = 0; ri<NUM_RELATIONS; ++ri) {
        for (uint32_t ci=0; ci<gSchema[ri]; ++ci) {
            processEqualityQueries(ri, ci);
        }
    }
    */
#ifdef LPDEBUG
    LPTimer.validationsProcessingIndex += LPTimer.getChrono(qproc);
#endif
}
void processPendingValidationsTask(uint32_t nThreads, uint32_t tid, void *args) {
    (void)tid; (void)nThreads; (void)args;// to avoid unused warning

    processEqualityQ(nThreads, tid, args);

    uint64_t totalPending = gPendingValidations.size();
    // get a validation ID - atomic operation
    for (uint64_t vi = gNextPending++; vi < totalPending; vi=gNextPending++) {
        auto& v = gPendingValidations[vi];
        uint64_t resPos = v.validationId - resIndexOffset;
        auto& atoRes = gPendingResults[resPos];
        if (atoRes) continue;
        if(isValidationConflict(v)) { atoRes = true; }
        //atoRes = isValidationConflict(v);
    } // while true take more validations 
}

