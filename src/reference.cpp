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

#include "include/DoubleIter.hpp"
#include "include/aligned_allocator.hpp"
#include "include/LPUtils.hpp"
#include "include/ReaderIO.hpp"
#include "include/LPThreadpool.hpp"
#include "include/SingleTaskPool.hpp"
#include "include/MultiTaskPool.hpp"

#include <iostream>
#include <fstream>
#include <ios>
#include <cstdio>
#include <iterator>
#include <map>
#include <unordered_map>
#include <set>
#include <vector>
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


#ifdef LPOMPTBB
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include "include/TBBUtils.hpp"

#include <omp.h>
#endif




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
using vector_a = std::vector<T, aligned_allocator<T, 16>>;
//using vector_a = std::vector<T>;

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

struct ColumnTransaction_t {
    vector_a<uint64_t> values;
    vector_a<tuple_t> tuples;
    uint64_t trans_id;
    ColumnTransaction_t(uint64_t tid) : values(vector_a<uint64_t>()), tuples(vector_a<tuple_t>()), trans_id(tid) {}
} ALIGNED_DATA;

typedef pair<uint64_t, tuple_t> Metadata_t;
struct ColumnStruct {
    // the trans_id the transactions are updated to inclusive

    vector<uint64_t> values;
    vector<Metadata_t> metadata;

    //vector<ColumnTransaction_t> transactions;
    //vector<uint64_t> transactionsORs;
    uint64_t transTo;


    
    char padding[8]; //for false sharing
    
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
/*
struct CTRSLessThan_t {
    ALWAYS_INLINE bool operator() (const ColumnTransaction_t& left, const ColumnTransaction_t& right) {
        return left.first < right.first;
    }
    ALWAYS_INLINE bool operator() (const ColumnTransaction_t& o, uint64_t target) {
        return o.first < target;
    }
    ALWAYS_INLINE bool operator() (uint64_t target, const ColumnTransaction_t& o) {
        return target < o.first;
    }
} CTRSLessThan;
*/
// transactions in each relation column - all tuples of same transaction in one vector

struct RelationColumns {
    std::unique_ptr<ColumnStruct[]> columns;
};
static std::unique_ptr<RelationColumns[]> gRelColumns;
static vector<uint32_t> gRequiredColumns;

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
// The general structure for each relation
struct RelationStruct {
    vector<pair<uint64_t, vector<tuple_t>>> transLogTuples;
    vector<unique_ptr<RelTransLog>> transLog;
    btree::btree_map<uint32_t, pair<uint64_t, uint64_t*>> insertedRows;
    
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
        else return l.isDelOp & !r.isDelOp;
    }
} TRMapPhaseByTrans;

typedef std::mutex RelTransLock;
//typedef LPSpinLock RelTransLock;
static unique_ptr<RelTransLock[]> gRelTransMutex;
static unique_ptr<vector<TRMapPhase>[]> gTransParseMapPhase;

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

/////////////////////////////////////////// STRUCTURES FOR STATS
/*
typedef pair<bool, uint32_t> SColType;
struct StatStruct {
    // columns info that appear as 1st predicates - bool=True means equality, False anything else
    vector<SColType> reqCols; 
};
struct StatCompEq_t {
    inline bool operator() (const SColType& l, const SColType& r) {
        return l.second == r.second;
    }
} StatCompEq;
struct StatComp_t {
    inline bool operator() (const SColType& l, const SColType& r) {
        return l.second < r.second;
    }
    inline bool operator() (const SColType& l, uint32_t target) {
        return l.second < target;
    }
    inline bool operator() (uint32_t target, const SColType& r) {
        return target < r.second;
    }
} StatComp;
struct StatComp2_t {
    inline bool operator() (const SColType& l, const SColType& r) {
        return lp::validation::unpackRel(l.second) < lp::validation::unpackRel(r.second);
    }
    inline bool operator() (const SColType& l, uint32_t target) {
        return lp::validation::unpackRel(l.second) < target;
    }
    inline bool operator() (uint32_t target, const SColType& r) {
        return target < lp::validation::unpackRel(r.second);
    }
} StatCompRel;
static unique_ptr<StatStruct[]> gStats;
//static vector<SColType> *gStatColumns;
*/

////////////////////////////////////////////////////

LPTimer_t LPTimer;

struct GlobalState {
    enum State : uint32_t { SCHEMA, TRANSACTION, VALIDATION, FORGET, FLUSH };
    State state;
    uint32_t nThreads;
} Globals;

//---------------------------------------------------------------------------

// JUST SOME FUNCTION DECLARATIONS THAT ARE DEFINED BELOW
static void checkPendingValidations(ISingleTaskPool*);
void processPendingValidationsTask(uint32_t nThreads, uint32_t tid, void *args);

static void processTransactionMessage(const Transaction& t, ReceivedMessage *msg); 
static inline void checkPendingTransactions(ISingleTaskPool *pool);
void processPendingIndexTask(uint32_t nThreads, uint32_t tid, void *args);

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
    //cerr << endl << "relations: " << NUM_RELATIONS << endl;
    const uint32_t rels = d.relationCount;
    for(uint32_t ri=0; ri<rels; ++ri) {
        //cerr << " " << gSchema[ci];
        gRelColumns[ri].columns.reset(new ColumnStruct[gSchema[ri]]);
        const uint32_t colsz = gSchema[ri];
        for (uint32_t ci=0; ci<colsz; ++ci)
            gRequiredColumns.push_back(lp::validation::packRelCol(ri, ci));
    }
}
//---------------------------------------------------------------------------

#ifdef LPDEBUG
static uint64_t gTotalTransactions = 0, gTotalTuples = 0, gTotalValidations = 0;
#endif
//---------------------------------------------------------------------------
static void processValidationQueries(const ValidationQueries& v, ReceivedMessage *msg, uint64_t tid = 0) {
#ifdef LPDEBUG
    auto start = LPTimer.getChrono();    
#endif 
    (void)tid;

    // TODO - OPTIMIZATION CAN BE DONE IF I JUST COPY THE WHOLE DATA instead of parsing it
    // try to put all the queries into a vector
    /*
    vector<LPQuery> queries;
    queries.reserve(v.queryCount);
    const char* qreader=v.queries;
    for (uint32_t i=0;i<v.queryCount;++i) {
        const Query *q=reinterpret_cast<const Query*>(qreader);
        queries.emplace_back(const_cast<Query*>(q));
        qreader+=sizeof(Query)+(sizeof(Query::Column)*q->columnCount);
    }
    */
    //cerr << v.validationId << "====" << v.from << ":" << v.to << "=" << v.queryCount << endl;
    //gPendingValidationsMutex.lock();
    //gPendingValidations.emplace_back(v.validationId, v.from, v.to, msg, move(queries));    
    //gPendingValidations.emplace_back(v.validationId, v.from, v.to, msg, move(vector<LPQuery>()));    
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


static atomic<uint64_t> gNextFRel;
void processForgetThreaded(uint32_t nThreads, uint32_t tid, void *args) {
    (void)tid; (void)nThreads; (void)args;// to avoid unused warning
    //cerr << "::: tid " << tid << "new" << endl;
    //auto& f = gF;
    auto f = *reinterpret_cast<Forget*>(args);
    uint64_t trans_id = f.transactionId;
    for (uint64_t ri = gNextFRel++; ri < NUM_RELATIONS; ri=gNextFRel++) {     
        auto& cRelCol = gRelColumns[ri];
        // clean the index columns
        /*
        const uint64_t colsz = gSchema[ri];
        for (uint32_t ci=0; ci<colsz; ++ci) {
            auto& cCol = cRelCol.columns[ci];
            auto ub = upper_bound(cCol.transactions.begin(), cCol.transactions.end(),
                        f.transactionId,
                        [](const uint64_t target, const ColumnTransaction_t& ct){ return target < ct.first; });
            
            cCol.transactions.erase(cCol.transactions.begin(), ub);
            cCol.transactionsORs.erase(cCol.transactionsORs.begin(), cCol.transactionsORs.begin()+(ub-cCol.transactions.begin()));
        }*/
        for (uint32_t ci=0; ci<gSchema[ri]; ++ci) {
            auto& cCol = cRelCol.columns[ci];
            if (cCol.values.empty()) continue;
            auto& colValues = cCol.values;
            auto& colMetadata = cCol.metadata;
            
            auto itbeg = SIter<Metadata_t, uint64_t>(colMetadata.data(), colValues.data());
            auto itend = SIter<Metadata_t, uint64_t>(colMetadata.data()+colMetadata.size(), colValues.data() + colValues.size());
            auto it = std::remove_if(itbeg, itend,
                            [=](SIter<Metadata_t, uint64_t>::reference meta) { return meta.a->first <= trans_id; });
            size_t delPos = it - itbeg;
            colMetadata.erase(colMetadata.begin()+delPos, colMetadata.end());
            colValues.erase(colValues.begin()+delPos, colValues.end());
        }
        // clean the transactions log 
        //auto& transLog = gRelations[ri].transLog; 
        //cerr << "size bef: " << transLog.size() << endl;
        //for (auto it = transLog.begin(), tend=transLog.end(); it!=tend && ((*it)->trans_id <= f.transactionId); ) {
        //    if ((*it)->aliveTuples == 0 && (*it)->last_del_id <= f.transactionId) { it = transLog.erase(it); tend=transLog.end(); }
        //    else ++it;
        //}
        
        // delete the transLogTuples
        auto& transLogTuples = gRelations[ri].transLogTuples;
        transLogTuples.erase(transLogTuples.begin(), 
                upper_bound(transLogTuples.begin(), transLogTuples.end(), f.transactionId,
                    [](const uint64_t target, const pair<uint64_t, vector<tuple_t>>& o){ return target < o.first; })
                );
    }
}

        
static void processForget(const Forget& f, ISingleTaskPool* pool) {
#ifdef LPDEBUG
    auto start = LPTimer.getChrono();
#endif
    
    gNextFRel = 0; 
    pool->startSingleAll(processForgetThreaded, (void*)&f);
    pool->waitSingleAll();
    
    (void)pool; (void)f;
/*
    // delete the transactions from the columns index
    for (uint32_t i=0; i<NUM_RELATIONS; ++i) {
        auto& cRelCol = gRelColumns[i];
        
        // clean the index columns
        for (uint32_t ci=0; ci<gSchema[i]; ++ci) {
            auto& cCol = cRelCol.columns[ci];
            auto& colValues = cCol.values;
            auto& colMetadata = cCol.metadata;
            std::sort(SIter<Metadata_t, uint64_t>(colMetadata.data(), colValues.data()), 
                SIter<Metadata_t, uint64_t>(colMetadata.data()+colMetadata.size(), colValues.data() + colValues.size()));
            
            auto ub = upper_bound(colMetadata.begin(), colMetadata.end(),
                        f.transactionId,
                        [](const uint64_t target, const Metadata_t& ct){ return target < ct.first; });
            
            colMetadata.erase(colMetadata.begin(), ub);
            colValues.erase(colValues.begin(), colValues.begin()+(ub-colMetadata.begin()));
        }

        // clean the transactions log
        //auto& transLog = gRelations[i].transLog;         
        //cerr << "size bef: " << transLog.size() << endl;
        //for (auto it = transLog.begin(), tend=transLog.end(); it!=tend && ((*it)->trans_id <= f.transactionId); ) {
        //    if ((*it)->aliveTuples == 0 && (*it)->last_del_id <= f.transactionId) { it = transLog.erase(it); tend=transLog.end(); }
        //    else ++it;
        //}

        // delete the transLogTuples
        auto& transLogTuples = gRelations[i].transLogTuples;
        transLogTuples.erase(transLogTuples.begin(), 
                upper_bound(transLogTuples.begin(), transLogTuples.end(), f.transactionId,
                    [](const uint64_t target, const pair<uint64_t, vector<tuple_t>>& o){ return target < o.first; })
                );
        //cerr << "size after: " << transLog.size() << endl;
    }
*/
#ifdef LPDEBUG
    LPTimer.forgets += LPTimer.getChrono(start);
#endif
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/////////////////// MAIN-READING STRUCTURES ///////////////////////


void inline parseValidation(uint32_t nThreads, uint32_t tid, void *args) {
    (void)tid; (void)nThreads;
#ifdef LPDEBUG
    auto start = LPTimer.getChrono();
#endif
    ReceivedMessage *msg = reinterpret_cast<ReceivedMessage*>(args);
    processValidationQueries(*reinterpret_cast<const ValidationQueries*>(msg->data.data()), msg, tid); 
    // TODO - do not delete the vector since it contains the data for the validation being used
    //delete msg;
#ifdef LPDEBUG
    LPTimer.validations += LPTimer.getChrono(start);
#endif
}




#ifdef LPOMPTBB
void inline initOpenMP(uint32_t nThreads) {
    //omp_set_dynamic(0);           // Explicitly disable dynamic teams
    omp_set_num_threads(nThreads);  // Use 4 threads for all consecutive parallel regions
}
void inline initTBB(uint32_t nThreads) {
   (void)nThreads;
   //tbb::task_scheduler_init init(tbb::task_scheduler_init::automatic); 
   //tbb::task_scheduler_init init(nThreads); 
}
#endif 

int main(int argc, char**argv) {
    uint64_t numOfThreads = 1;
    if (argc > 1) {
        numOfThreads = strtol(argv[1], NULL, 10);
        if (numOfThreads == LONG_MAX)
            numOfThreads = 1;
    }
    //cerr << "Number of threads: " << numOfThreads << endl;
    Globals.nThreads = numOfThreads;

#ifdef LPOMPTBB
    initOpenMP(numOfThreads);
    initTBB(numOfThreads);
#endif

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
    }

    // do some initial reserves or initializations
    gPendingValidations.reserve(2048); 
    for (uint32_t i=0; i<NUM_RELATIONS; ++i) gTransParseMapPhase[i].reserve(512);
    for (uint32_t i=0; i<NUM_RELATIONS; ++i) {
        gRelations[i].transLog.reserve(512);
        gRelations[i].transLogTuples.reserve(1024);
    }
    // allocate global structures based on thread number
    //gStats.reset(new StatStruct[numOfThreads+1]);

    // allocate the workers
    SingleTaskPool workerThreads(numOfThreads, processPendingValidationsTask);
    //SingleTaskPool workerThreads(1, processPendingValidationsTask);
    workerThreads.initThreads();
    // leave two available workes - master - Reader
    //MultiTaskPool multiPool(std::max(numOfThreads-4, (uint64_t)2));
    //MultiTaskPool multiPool(1);
    //multiPool.initThreads();
    //multiPool.startAll();

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
                        //multiPool.addTask(parseValidation, static_cast<void*>(msg)); 
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
                    checkPendingValidations(&workerThreads);
                    Globals.state = GlobalState::FLUSH;
                    processFlush(*reinterpret_cast<const Flush*>(msg->data.data()), isTestdriver); 
                    delete msg;
                    break;
                case MessageHead::Forget: 
#ifdef LPDEBUG
                    ++totalForgets; 
#endif
                    // check if we have pending transactions to be processed
                    checkPendingValidations(&workerThreads);
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
                        //multiPool.destroy();
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
            //for (uint32_t c=0; c<o.rowCount; ++c) ptr[c] = o.keys[c];
            memcpy(ptr, o.keys, sizeof(uint64_t)*o.rowCount);
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
            //uint32_t sz = relCols*o.rowCount;
            //for (uint32_t c=0; c<sz; ++c) tptr[c] = o.values[c];
            memcpy(tptr, o.values, sizeof(uint64_t)*relCols*o.rowCount);
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


/*
   struct TRMapPhase {
   uint64_t trans_id;
   bool isDelOp;
   uint32_t rowCount;
   uint64_t *values; // delete op => row keys to delete | insert => tuples
   }
//static unique_ptr<vector<TRMapPhase>[]> gTransParseMapPhase;
 */
/*
//static void updateRequiredColumns(uint64_t ri, vector<SColType>::iterator colBegin, vector<SColType>::iterator colEnd) {
static void updateRequiredColumns(uint64_t ri) {
    // PHASE TWO OF THE ALGORITHM IN THIS STAGE IS TO INCREMENTALLY UPDATE
    // THE INDEXES ONLY FOR THE COLUMNS THAT ARE GOING TO BE REQUESTED IN THE 
    // FOLOWING VALIDATION SESSION - 1st predicates only for now
    //for (SColType& cp : *statCols) cerr << "is Op::Equal " << cp.first << " col: " << cp.second << endl; 
    auto& relation = gRelations[ri];
    auto& relColumns = gRelColumns[ri].columns;
    
    uint64_t updatedUntil = relColumns[0].transTo;
    // Use lower_bound to automatically jump to the transaction to start
    auto transFrom = lower_bound(relation.transLogTuples.begin(), relation.transLogTuples.end(), updatedUntil, TransLogComp);
    auto tEnd=relation.transLogTuples.end();
    // for each column to be indexed
    uint32_t colsz=gSchema[ri];
//#pragma omp parallel for schedule(static, 1) num_threads(4)
    for (uint32_t col=0; col<colsz; ++col) {
        //tbb::parallel_for ((uint32_t)0, gSchema[ri], [&] (uint32_t col) {
    //tbb::parallel_for (tbb::blocked_range<uint32_t>(0, gSchema[ri], 20), [&] (const tbb::blocked_range<uint32_t>& r) {
    //    for (uint32_t col=r.begin(); col<r.end(); ++col) {

        //uint32_t rel,col;
        //for (; colBegin!=colEnd; ++colBegin) {
        //    lp::validation::unpackRelCol(colBegin->second, rel, col);
        //cerr << "relation: " << ri << " got rel " << rel << " col " << col << endl;
        auto& colTransactions = relColumns[col].transactions;
        auto& colTransactionsORs = relColumns[col].transactionsORs;

        // for all the transactions in the relation
        for(auto trp=transFrom; trp!=tEnd; ++trp) {
            colTransactionsORs.push_back(0);
            // allocate vectors for the current new transaction to put its data
            colTransactions.emplace_back(trp->first, move(vector<CTransStruct>()));
            colTransactions.back().second.reserve(trp->second.size());
            auto& vecBack = colTransactions.back().second;
            for (auto tpl : trp->second) {
                vecBack.emplace_back(tpl[col], tpl);
                colTransactionsORs.back() |= tpl[col];
            }
            sort(vecBack.begin(), vecBack.end(), ColTransValueLess);
            //cerr << "OR: " << colTransactionsORs.back() << endl;
            // add the sentinel value
            //vecBack.emplace_back(UINT64_MAX, nullptr);
        }
        if(!relation.transLogTuples.empty())
            relColumns[col].transTo = max(relation.transLogTuples.back().first+1, updatedUntil);
        //cerr << "col " << col << " ends to " << relColumns[col].transTo << endl;
      //      }});
    }
}
*/
static std::atomic<uint64_t> gNextIndex;

void processPendingIndexTask(uint32_t nThreads, uint32_t tid, void *args) {
    (void)tid; (void)nThreads; (void)args;// to avoid unused warning
    //cerr << "::: tid " << tid << "new" << endl;

    for (uint64_t ri = gNextIndex++; ri < NUM_RELATIONS; ri=gNextIndex++) {
        //#pragma omp parallel for schedule(static, 1)  
        //for(uint64_t ri = 0; ri < NUM_RELATIONS; ++ri) {

        //auto colpair = std::equal_range(gStatColumns->begin(), gStatColumns->end(), ri, StatCompRel);
        //auto colBegin = colpair.first, colEnd = colpair.second; 

        // take the vector with the transactions and sort it by transaction id in order to apply them in order
        auto& relTrans = gTransParseMapPhase[ri];
        if (unlikely(relTrans.empty())) { 
            // TODO - we have to run this regardless of transactions since some
            // columns might have to use previous transactions and be called for the first time
            //updateRequiredColumns(ri, colBegin, colEnd);
            //updateRequiredColumns(ri);
            continue; 
        }

        //cerr << "tid " << tid << " got " << ri << " = " << relTrans.size() << endl;

        //std::sort(relTrans.begin(), relTrans.end(), TRMapPhaseByTrans);

        auto& relation = gRelations[ri];
        //auto& relColumns = gRelColumns[ri].columns;
        uint32_t relCols = gSchema[ri];

        uint64_t lastTransId = relTrans[0].trans_id;

        // for each transaction regarding this relation
        vector<tuple_t> operations;
        operations.reserve(512);
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

        // update with new transactions
        //updateRequiredColumns(ri);
    } // end of while true
}

static std::atomic<uint32_t> gNextReqCol;

void ALWAYS_INLINE updateIndexRelCol(uint32_t tid, uint32_t ri, uint32_t col) { (void)tid;
    auto& relation = gRelations[ri];
    auto& relColumns = gRelColumns[ri].columns;
    
    uint64_t updatedUntil = relColumns[col].transTo;
    if (relation.transLogTuples.empty() || lp_EQUAL(updatedUntil, relation.transLogTuples.back().first)) return;
    
    // Use lower_bound to automatically jump to the transaction to start
    auto transFrom = lower_bound(relation.transLogTuples.begin(), relation.transLogTuples.end(), updatedUntil, TransLogComp);
    auto tEnd=relation.transLogTuples.end();
    //auto& colTransactions = relColumns[col].transactions;
    //auto& colTransactionsORs = relColumns[col].transactionsORs;

    auto& colValues = relColumns[col].values;
    auto& colMetadata = relColumns[col].metadata;
    
    // for all the transactions in the relation
    for(auto trp=transFrom; trp!=tEnd; ++trp) {
        // allocate vectors for the current new transaction to put its data
        //colTransactions.emplace_back(trp->first);
        //auto& values = colTransactions.back().values;
        //auto& tuples = colTransactions.back().tuples;
        //const unsigned int trpsz = trp->second.size();
        //values.reserve(trpsz);
        //tuples.reserve(trpsz);
        //colTransactionsORs.push_back(0);
        for (auto tpl : trp->second) {
            //values.push_back(tpl[col]);
            //tuples.push_back(tpl);
            //colTransactionsORs.back() |= tpl[col];
            colValues.push_back(tpl[col]);
            colMetadata.emplace_back(trp->first, tpl);
        }

        //std::sort(SIter<uint64_t, tuple_t>(values.data(), tuples.data()), 
        //      SIter<uint64_t, tuple_t>(values.data()+trpsz, tuples.data()+trpsz));
    }
    // no need to check for empty since now we update all the columns and there is a check for emptyness above
    relColumns[col].transTo = max(relation.transLogTuples.back().first+1, updatedUntil);

    // TODO - MANY MANY - TODO - MANY THINGS TO DO HERE - Btree - OR UPDATE INCRMENETALLY
    // TODO - OR USER INPLACE_MERGE of std::
    std::sort(SIter<uint64_t, pair<uint64_t, tuple_t>>(colValues.data(), colMetadata.data()), 
          SIter<uint64_t, pair<uint64_t, tuple_t>>(colValues.data()+colValues.size(), colMetadata.data()+colMetadata.size()));
}

void processUpdateIndexTask(uint32_t nThreads, uint32_t tid, void *args) {
    (void)tid; (void)nThreads; (void)args;// to avoid unused warning
    uint64_t totalCols = gRequiredColumns.size();
    for (uint64_t rc = gNextReqCol++; rc < totalCols; rc=gNextReqCol++) {
        uint32_t ri, col;
        lp::validation::unpackRelCol(gRequiredColumns[rc], ri, col);
        updateIndexRelCol(tid, ri, col);
    } // end of while columns to update     
}

static inline void checkPendingTransactions(ISingleTaskPool *pool) {
#ifdef LPDEBUG
    auto startIndex = LPTimer.getChrono();
#endif
    //cerr << "::: session start ::::" << endl;
    gNextIndex = 0;
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
    //processUpdateIndexTask(0, 0, nullptr);
    pool->startSingleAll(processUpdateIndexTask);
    pool->waitSingleAll();
#ifdef LPDEBUG
    LPTimer.updateIndex += LPTimer.getChrono(startUpdIndex);
#endif
}


static uint64_t resIndexOffset = 0;
static std::atomic<uint64_t> gNextPending;

static void checkPendingValidations(ISingleTaskPool *pool) {
    if (unlikely(gPendingValidations.empty())) return;

    // check if there is any pending index creation to be made before checking validation
    checkPendingTransactions(pool);

#ifdef LPDEBUG
    auto start = LPTimer.getChrono();
#endif

    // find the MIN validation ID to coordinate the indexing of the results
    //resIndexOffset = UINT64_MAX;
    //for (auto& pv : gPendingValidations) if (pv.validationId < resIndexOffset) resIndexOffset = pv.validationId;
    resIndexOffset = gPendingValidations[0].validationId; // only master handles the messages now so they are in order
    
    auto gPRsz = gPendingResults.size();
    if (gPVunique > gPRsz)
        gPendingResults.resize(gPVunique);
    memset(gPendingResults.data(), 0, sizeof(PendingResultType)*gPRsz);
    //std::fill(gPendingResults.begin(), gPendingResults.end(), 0);
    gNextPending = 0;

    // sort the validations by query count in order to start the heavy ones earlier
    //std::sort(gPendingValidations.begin(), gPendingValidations.end(), 
      //      [](const LPValidation& left, const LPValidation& right){ return left.queryCount > right.queryCount; });
    //std::sort(gPendingValidations.begin(), gPendingValidations.end(), LPValCompQCount);

    pool->startSingleAll(processPendingValidationsTask);
    pool->waitSingleAll();

    //(void)pool;
    //processPendingValidationsTask(0, 0);

    // update the results - you can get the validation id by adding resIndexOffset to the position
    for (uint64_t i=0, valId=resIndexOffset; i<gPVunique; ++i, ++valId) { 
        gQueryResults.emplace_back(valId, gPendingResults[i]);
    }
    gPendingValidations.clear();
    gPVunique = 0;
#ifdef LPDEBUG
    LPTimer.validationsProcessing += LPTimer.getChrono(start);
#endif
}

//////////////////////////////////////////////////////////
// VALIDATION
//////////////////////////////////////////////////////////

//typedef CTransStruct TupleType;
typedef tuple_t TupleType;
//typedef vector<Query::Column>::iterator PredIter;
typedef Query::Column* PredIter;

bool ALWAYS_INLINE isTupleConflict(PredIter cbegin, PredIter cend, const TupleType& tup) {
    for (; cbegin<cend; ++cbegin) {
        auto& c = *cbegin;
        // make the actual check
        switch (c.op) {
            case Op::Equal: 
                if(!lp_EQUAL(tup[c.column], c.value)) return false; 
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
                if(lp_EQUAL(tup[c.column], c.value)) return false; 
                break;
        } 
    } // end of single query predicates
    return true;    
}

//bool ALWAYS_INLINE isTupleRangeConflict(vector_a<TupleType>::const_iterator tupFrom, vector_a<TupleType>::const_iterator tupTo, PredIter cbegin, PredIter cend) {
bool ALWAYS_INLINE isTupleRangeConflict(TupleType *tupFrom, TupleType *tupTo, PredIter cbegin, PredIter cend) {
    for(; tupFrom!=tupTo; ++tupFrom) {  
        if (isTupleConflict(cbegin, cend, *tupFrom)) return true;
    } // end of all tuples for this transaction
    return false;
    //return std::find_if(tupFrom, tupTo, [=](TupleType& tup) { return isTupleConflict(cbegin, cend, tup);}) != tupTo;
}

//static bool inline isTransactionConflict(vector_a<CTransStruct>& transValues, Column pFirst) {
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

auto kernelZero = [](uint64_t t) { return t != 0; };

//vector<uint64_t> cres;
//vector<uint64_t> resTuples;
/*
bool isTupleRangeConflict(TupleType *tupFrom, TupleType *tupTo, 
        PredIter cbegin, PredIter cend, ColumnStruct *relColumns, unsigned int pos) {
    // copy the eligible tuples into our mask vector
    //vector<uint64_t> cres; 
    size_t csz;
    vector<uint64_t> resTuples; 
    //resTuples.resize(0);
    resTuples.reserve(tupTo-tupFrom);
    for (; tupFrom!=tupTo; ++tupFrom) resTuples.push_back((uint64_t)*tupFrom);
    size_t activeSize = resTuples.size();
    for (; cbegin<cend; ++cbegin) {
        //cerr << "active: " << activeSize << endl;
        auto& c = *cbegin;    
        auto& cTransactions = relColumns[c.column].transactions[pos];
        auto& transValues = cTransactions.values;
        auto& transTuples = cTransactions.tuples;
        decltype(transValues.begin()) tBegin = transValues.begin(), tEnd=transValues.end();
        size_t tupFromIdx{0}, tupToIdx{transTuples.size()};
        switch (c.op) {
            case Op::Equal: 
                {
                    if (transValues[0] > c.value || transValues.back() < c.value) return false;
                    auto tp = std::equal_range(tBegin, tEnd, c.value);
                    if (tp.second == tp.first) return false;
                    tupFromIdx = (tp.first - tBegin); tupToIdx = tupFromIdx + (tp.second-tp.first);
                    break;}
            case Op::Less: 
                if (transValues[0] >= c.value) return false;
                tupToIdx -= (tEnd-std::lower_bound(tBegin, tEnd, c.value));
                if (tupToIdx == tupFromIdx) return false;
                break;
            case Op::LessOrEqual: 
                if (transValues[0] > c.value) return false;
                tupToIdx -= (tEnd-std::upper_bound(tBegin, tEnd, c.value)); 
                if (tupToIdx == tupFromIdx) return false;
                break;
            case Op::Greater: 
                if (transValues.back() <= c.value) return false;
                tupFromIdx += (std::upper_bound(tBegin, tEnd, c.value)-tBegin);
                if (tupToIdx == tupFromIdx) return false;
                break;
            case Op::GreaterOrEqual: 
                if (transValues.back() < c.value) return false;
                tupFromIdx += (std::lower_bound(tBegin, tEnd, c.value)-tBegin);
                if (tupToIdx == tupFromIdx) return false;
                break;
            default: 
                // check if the active tuples have a value != to the predicate
                for (uint64_t *resPtr=resTuples.data(), *resend = resTuples.data()+activeSize; resPtr<resend;) {
                    if (((tuple_t)*resPtr++)[c.column] == c.value) *(resPtr-1) = 0;
                }
                goto LBL_CHECK_END;
        }

        // this check is done for all the operators apart from !=
        // we have to check if the active tuples are inside the result set returned
        csz = tupToIdx-tupFromIdx;
        if (csz > 512 && activeSize > 512) {
            //cerr << "csz: " << csz << " active: " << activeSize << endl;
            //cres.resize(0);
            //for (size_t i=tupFromIdx; i<tupToIdx; ++i) cres.push_back((uint64_t)transTuples[i]);
            //cres.resize(csz);
            //memcpy(cres.data(), transTuples.data()+tupFromIdx, sizeof(tuple_t)*csz);
            //std::sort(cres.begin(), cres.end());
            uint64_t *cres = (uint64_t*)alloca(csz*sizeof(uint64_t));
            memcpy(cres, transTuples.data()+tupFromIdx, sizeof(tuple_t)*csz);
            std::sort(cres, cres+csz);
            const size_t extra = activeSize & 1;
            uint64_t *resPtr = resTuples.data();
            const uint64_t *cresb = cres;
            const uint64_t *crese = cres + csz;
            //for (size_t i=0, nsz=activeSize-extra; i<nsz; i+=2) {
            for (auto nsz=resPtr+activeSize-extra; resPtr<nsz; resPtr += 2) {
                //if (!lp::utils::binary_cmov(cresb, csz, *resPtr)) *resPtr = 0;
                //if (!lp::utils::binary_cmov(cresb, csz, *(resPtr + 1))) *(resPtr+1) = 0; 
                if (!std::binary_search(cresb, crese, *resPtr)) *resPtr = 0;
                if (!std::binary_search(cresb, crese, *(resPtr + 1))) *(resPtr+1) = 0; 
            }
            //if (extra && !lp::utils::binary_cmov(cresb, csz, *resPtr)) *resPtr = 0;
            if (extra && !std::binary_search(cresb, crese, *resPtr)) *resPtr = 0;
        } else {
            const uint64_t *transPtr = (uint64_t*)(transTuples.data()+tupFromIdx);
            const size_t extra = activeSize & 1;
            uint64_t *resPtr = resTuples.data();
            const uint64_t *resEnd = resPtr+activeSize-extra;
            //for (size_t i=0; i<activeSize; ++i) {
            for (; resPtr<resEnd; resPtr += 2) {
                //auto start = LPTimer.getChrono();
                if (!lp::simd::exists_avx(transPtr, csz, *resPtr)) *resPtr = 0;
                if (!lp::simd::exists_avx(transPtr, csz, *(resPtr+1))) *(resPtr+1) = 0; 
                //gTimeSearch += LPTimer.getChrono(start);
            }
            if (extra && !lp::simd::exists_avx(transPtr, csz, *resPtr)) *resPtr = 0;
        }
LBL_CHECK_END:
        // check if we have any valid tuple left otherwise return false
        activeSize = std::partition(resTuples.data(), resTuples.data()+activeSize, kernelZero) - resTuples.data();
        //cerr << "active (after): " << activeSize << endl;
        //if (activeSize == 0) return false;
        if (activeSize < 172) return isTupleRangeConflict(reinterpret_cast<tuple_t*>(resTuples.data()), reinterpret_cast<tuple_t*>(resTuples.data()+activeSize), ++cbegin, cend);
    }
    return true;
}

static bool isTransactionConflict(const ColumnTransaction_t& transaction, Column pFirst, PredIter cbegin, PredIter cend, ColumnStruct *relColumns, unsigned int pos) {
    auto& transValues = transaction.values;
    auto& transTuples = transaction.tuples;
    decltype(transValues.begin()) tBegin = transValues.begin(), tEnd=transValues.end();
    TupleType *tupFrom{const_cast<tuple_t*>(transTuples.data())}, *tupTo{const_cast<tuple_t*>(transTuples.data()+transTuples.size())};
    // find the valid tuples using range binary searches based on the first predicate
    switch (pFirst.op) {
        case Op::Equal: 
            {
                if (transValues[0] > pFirst.value || transValues.back() < pFirst.value) return false;
                auto tp = std::equal_range(tBegin, tEnd, pFirst.value);
                if (tp.second == tp.first) return false;
                tupFrom += (tp.first - tBegin); tupTo -= (tEnd-tp.second);
                break;}
        case Op::Less: 
            if (transValues[0] >= pFirst.value) return false;
            tupTo -= (tEnd-std::lower_bound(tBegin, tEnd, pFirst.value));
            if (tupTo == tupFrom) return false;
            break;
        case Op::LessOrEqual: 
            if (transValues[0] > pFirst.value) return false;
            tupTo -= (tEnd-std::upper_bound(tBegin, tEnd, pFirst.value)); 
            if (tupTo == tupFrom) return false;
            break;
        case Op::Greater: 
            if (transValues.back() <= pFirst.value) return false;
            tupFrom += (std::upper_bound(tBegin, tEnd, pFirst.value)-tBegin);
            if (tupTo == tupFrom) return false;
            break;
        case Op::GreaterOrEqual: 
            if (transValues.back() < pFirst.value) return false;
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
    
    //if (tupTo - tupFrom < 128) return isTupleRangeConflict(tupFrom, tupTo, cbegin, cend);
    if (tupTo - tupFrom < 64) return isTupleRangeConflict(tupFrom, tupTo, cbegin, cend);
    else return isTupleRangeConflict(tupFrom, tupTo, cbegin, cend, relColumns, pos);
    //return isTupleRangeConflict(tupFrom, tupTo, cbegin, cend, relColumns, pos);
}
*/

static bool isConflict(LPValidation& v, Column pFirst, PredIter cbegin, PredIter cend, ColumnStruct *relColumns) {
    auto& colValues = relColumns[pFirst.column].values;
    auto& colMetadata = relColumns[pFirst.column].metadata;
    decltype(colValues.begin()) tBegin = colValues.begin(), tEnd=colValues.end();
    Metadata_t *tupFrom{const_cast<Metadata_t*>(colMetadata.data())}, *tupTo{const_cast<Metadata_t*>(colMetadata.data()+colMetadata.size())};
    // find the valid tuples using range binary searches based on the first predicate
    switch (pFirst.op) {
        case Op::Equal: 
            {
                if (colValues[0] > pFirst.value || colValues.back() < pFirst.value) return false;
                auto tp = std::equal_range(tBegin, tEnd, pFirst.value);
                if (tp.second == tp.first) return false;
                tupFrom += (tp.first - tBegin); tupTo -= (tEnd-tp.second);
                break;}
        case Op::Less: 
            if (colValues[0] >= pFirst.value) return false;
            tupTo -= (tEnd-std::lower_bound(tBegin, tEnd, pFirst.value));
            if (tupTo == tupFrom) return false;
            break;
        case Op::LessOrEqual: 
            if (colValues[0] > pFirst.value) return false;
            tupTo -= (tEnd-std::upper_bound(tBegin, tEnd, pFirst.value)); 
            if (tupTo == tupFrom) return false;
            break;
        case Op::Greater: 
            if (colValues.back() <= pFirst.value) return false;
            tupFrom += (std::upper_bound(tBegin, tEnd, pFirst.value)-tBegin);
            if (tupTo == tupFrom) return false;
            break;
        case Op::GreaterOrEqual: 
            if (colValues.back() < pFirst.value) return false;
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
    
    for(; tupFrom!=tupTo; ++tupFrom) {  
        if (tupFrom->first>=v.from && tupFrom->first<=v.to && isTupleConflict(cbegin, cend, tupFrom->second)) return true;
    } // end of all tuples for this transaction
    return false;
}


static bool isValidationConflict(LPValidation& v) {
    // TODO - MAKE A PROCESSING OF THE QUERIES AND PRUNE SOME OF THEM OUT
    const ValidationQueries& vq = *reinterpret_cast<ValidationQueries*>(v.rawMsg->data.data());
    /*
    cerr << "\n========= validation " << v.validationId << " =========" << endl; 
    cerr << "qc: " << vq.queryCount << " from: " << vq.from << " to: " << vq.to << endl; 
    */
    const char* qreader = vq.queries;
    uint32_t columnCount;
   /*
    vector<pair<uint32_t, Column>> allPreds;
    for (uint32_t i=0; i<vq.queryCount; ++i, qreader+=sizeof(Query)+(sizeof(Query::Column)*columnCount)) {
        const Query *rq = (reinterpret_cast<const Query*>(qreader));
        columnCount = rq->columnCount;
        cerr << "-- query: rel: " << rq->relationId << " columns: " << rq->columnCount << endl;
        for (uint32_t j=0; j<rq->columnCount; ++j) {
            cerr << rq->columns[j];
            allPreds.emplace_back(i, rq->columns[j]);
        }
        cerr << endl;
    }
    qreader = vq.queries;

    auto comp = [](const pair<uint32_t, Column>& l, pair<uint32_t, Column>& r) {
                return ColumnCompCol(l.second, r.second);
            };
    auto compeq = [](const pair<uint32_t, Column>& l, pair<uint32_t, Column>& r) {
                return !ColumnCompCol(l.second, r.second) && !ColumnCompCol(l.second, r.second);
            };
    sort(allPreds.begin(), allPreds.end(), comp);
    auto uniqit = std::unique(allPreds.begin(), allPreds.end(), compeq);
    cerr << "uniques: " << std::distance(allPreds.begin(), uniqit) << " : " << allPreds.size() << endl;
    cerr << "\npredicates sorted :::" << endl;
    for (auto c : allPreds) {
        cerr << c.first << "_" << c.second << " ";
    }
    cerr << endl;

    cerr << "checked: " << endl;
    */
    for (uint32_t i=0; i<vq.queryCount; ++i, qreader+=sizeof(Query)+(sizeof(Query::Column)*columnCount)) {
        Query& rq=*const_cast<Query*>(reinterpret_cast<const Query*>(qreader));
        columnCount = rq.columnCount;
       
        //cerr << " " << i;

        if (unlikely(rq.columnCount == 0)) { 
            //cerr << "empty: " << v.validationId << endl; 
            auto& transactionsCheck = gRelations[rq.relationId].transLogTuples;
            auto transFromCheck = std::lower_bound(transactionsCheck.begin(), transactionsCheck.end(), v.from, TransLogComp);
            auto transToCheck = std::upper_bound(transFromCheck, transactionsCheck.end(), v.to, TransLogComp);
            if (transFromCheck == transToCheck) { 
                // no transactions exist for this query
                continue;
            } else { 
                // transactions exist for this query so it is a conflict
                return true;
            }; 
        }
        
        uint32_t colCountUniq = lp::query::preprocess(rq); 
        if (!lp::query::satisfiable(&rq, colCountUniq)) { /*cerr << "rej" << endl;*/ continue; } // go to the next query
        //cerr << "passed" << endl;
        
        auto cbegin = reinterpret_cast<Query::Column*>(rq.columns),
            cend = cbegin + colCountUniq;
        auto pFirst = *reinterpret_cast<Query::Column*>(rq.columns);
        // just find the range of transactions we want in this relation
        auto& relColumns = gRelColumns[rq.relationId].columns;
        //////////////////////////////

        auto& colValues = relColumns[pFirst.column].values;
        auto& colMetadata = relColumns[pFirst.column].metadata;
        if (colValues.empty()) continue; //next query since this one is not satisfiable

        if (isConflict(v, pFirst, cbegin+1, cend, relColumns.get())) return true;
        else continue;

        /////////////////////////////
        /*
        auto& transactions = relColumns[pFirst.column].transactions;
        auto transFrom = std::lower_bound(transactions.begin(), transactions.end(), v.from, CTRSLessThan);
        auto transTo = std::upper_bound(transFrom, transactions.end(), v.to, CTRSLessThan);
        
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
                if (    !((!(cb.op) & !lp_EQUAL((relColumns[cb.column].transactionsORs[pos] & cb.value), cb.value))
                        || (!(cb1.op) & !lp_EQUAL((relColumns[cb1.column].transactionsORs[pos] & cb1.value), cb1.value)))
                && isTransactionConflict(*transFrom, pFirst, cbSecond, cend, relColumns.get(), pos)) { return true; }
            } // end of all the transactions for this relation for this specific query
        } else {
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
        /*
        //cerr << ":: val " << v.validationId << endl;
        for(; transFrom<transTo; ++transFrom) {  
            if (isTransactionConflict(transFrom->second, pFirst)) { return true; }
        } // end of all the transactions for this relation for this specific query
        */
    }// end for all queries
    return false;
}

void processPendingValidationsTask(uint32_t nThreads, uint32_t tid, void *args) {
    (void)tid; (void)nThreads; (void)args;// to avoid unused warning

    uint64_t totalPending = gPendingValidations.size();
    // get a validation ID - atomic operation
    for (uint64_t vi = gNextPending++; vi < totalPending; vi=gNextPending++) {
    //for (uint64_t vi = gNextPending--; vi < totalPending; vi=gNextPending--) {
        //d:cerr << vi << " ";
        auto& v = gPendingValidations[vi];
        uint64_t resPos = v.validationId - resIndexOffset;
        auto& atoRes = gPendingResults[resPos];
        //if(isValidationConflict(v)) { atoRes = true; }
        atoRes = isValidationConflict(v);
        delete v.rawMsg;
    } // while true take more validations 
}

