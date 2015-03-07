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
#include "include/LPUtils.hpp"
#include "include/TBBUtils.hpp"
#include "include/ReaderIO.hpp"

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
#include "include/SingleTaskPool.hpp"
#include "include/MultiTaskPool.hpp"
#include "include/LPSpinLock.hpp"

#include "include/cpp_btree/btree_map.h"

#include "include/ReferenceTypes.hpp"
#include "include/LPQueryTypes.hpp"


#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include <omp.h>


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

// Custom data structures to hold data
struct CTransStruct {
    uint64_t value;
    tuple_t tuple;
    CTransStruct (uint64_t v, tuple_t t) : value(v), tuple(t) {}
};
struct CTRSValueLessThan_t {
    inline bool operator()(const CTransStruct& l, const CTransStruct& r) {
        return l.value < r.value;
    }
    inline bool operator() (const CTransStruct& o, uint64_t target) {
        return o.value < target;
    }
    inline bool operator() (uint64_t target, const CTransStruct& o) {
        return target < o.value;
    }
} ColTransValueLess;
std::ostream& operator<< (std::ostream& os, const CTransStruct& o) {
    os << "{" << "-" << o.value << "}";
    return os;
}
typedef pair<uint64_t, vector<CTransStruct>> ColumnTransaction_t;
struct ColumnStruct {
    // the trans_id the transactions are updated to inclusive
    uint64_t transTo;
    vector<ColumnTransaction_t> transactions;
    ColumnStruct() : transTo(0) {}
};
struct CTRSLessThan_t {
    inline bool operator() (const ColumnTransaction_t& left, const ColumnTransaction_t& right) {
        return left.first < right.first;
    }
    inline bool operator() (const ColumnTransaction_t& o, uint64_t target) {
        return o.first < target;
    }
    inline bool operator() (uint64_t target, const ColumnTransaction_t& o) {
        return target < o.first;
    }
} CTRSLessThan;

// transactions in each relation column - all tuples of same transaction in one vector

struct RelationColumns {
    std::unique_ptr<ColumnStruct[]> columns;
};
static std::unique_ptr<RelationColumns[]> gRelColumns;


struct RelTransLog {
    uint64_t trans_id;
    uint64_t last_del_id;
    uint64_t aliveTuples;
    uint64_t rowCount;
    //uint64_t *tuples;
    std::unique_ptr<uint64_t[]> tuples;

    RelTransLog(uint64_t tid, uint64_t* tpl, uint64_t _rowCount) : trans_id(tid), last_del_id(tid), aliveTuples(_rowCount), rowCount(_rowCount) {
        if (tpl != nullptr) tuples.reset(tpl);
    } 
};
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
};

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
};
struct TransMapPhase_t {
    inline bool operator()(const TRMapPhase& l, const TRMapPhase& r) {
        if (l.trans_id < r.trans_id) return true;
        else if (r.trans_id < l.trans_id) return false;
        else return l.isDelOp && !r.isDelOp;
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
static vector<SColType> *gStatColumns;

////////////////////////////////////////////////////

LPTimer_t LPTimer;

struct GlobalState {
    enum State : uint32_t { SCHEMA, TRANSACTION, VALIDATION, FORGET, FLUSH };
    State state;
    uint32_t nThreads;
} Globals;

//---------------------------------------------------------------------------

// JUST SOME FUNCTION DECLARATIONS THAT ARE DEFINED BELOW
static void checkPendingValidations(SingleTaskPool&);
static void processPendingValidationsTask(uint32_t nThreads, uint32_t tid);

static void processTransactionMessage(const Transaction& t, ReceivedMessage *msg);
static inline void checkPendingTransactions(SingleTaskPool& pool);
static void processPendingIndexTask(uint32_t nThreads, uint32_t tid);

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
    for(uint32_t ci=0; ci<d.relationCount; ++ci) {
        gRelColumns[ci].columns.reset(new ColumnStruct[gSchema[ci]]);
    }
}
//---------------------------------------------------------------------------

#ifdef LPDEBUG
static uint64_t gTotalTransactions = 0, gTotalTuples = 0, gTotalValidations = 0;
#endif
//---------------------------------------------------------------------------
static void processValidationQueries(const ValidationQueries& v, ReceivedMessage *msg, uint64_t tid = 0) {

#ifdef LPDEBUG
    //auto start = LPTimer.getChrono();    
#endif 
    (void)tid;

    // TODO - OPTIMIZATION CAN BE DONE IF I JUST COPY THE WHOLE DATA instead of parsing it
    // try to put all the queries into a vector
    vector<LPQuery> queries;
    queries.reserve(v.queryCount);
    const char* qreader=v.queries;
    for (uint32_t i=0;i<v.queryCount;++i) {
        const Query *q=reinterpret_cast<const Query*>(qreader);
        /*
           LPQuery nQ;
        //cerr << v.validationId << "====" << v.from << ":" << v.to << nQ << endl;
        if (likely(lp::query::parse(q, gSchema[q->relationId], &nQ))) {
        // this is a valid query
        nQ.relationId = q->relationId;

        //if (likely(!nQ.predicates.empty())) {
        //    // gather statistics    
        //    auto& p = nQ.predicates[0];
        //    uint32_t rc = lp::validation::packRelCol(nQ.relationId, p.column);
        //    gStats[tid].reqCols.emplace_back((p.op == Op::Equal), rc);
        //}

        queries.push_back(move(nQ));
        }
         */
        queries.emplace_back(const_cast<Query*>(q));
        qreader+=sizeof(Query)+(sizeof(Query::Column)*q->columnCount);
    }
    //  cerr << v.validationId << "====" << v.from << ":" << v.to << "=" << v.queryCount << "=" << queries << endl;
    gPendingValidationsMutex.lock();
    gPendingValidations.emplace_back(v.validationId, v.from, v.to, msg, move(queries));    
    // update the global pending validations to reflect this new one
    ++gPVunique;
    gPendingValidationsMutex.unlock();

#ifdef LPDEBUG
    //LPTimer.validations += LPTimer.getChrono(start);
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

static void processForget(const Forget& f) {
#ifdef LPDEBUG
    auto start = LPTimer.getChrono();
#endif

    // delete the transactions from the columns index
    for (uint32_t i=0; i<NUM_RELATIONS; ++i) {
        auto& cRelCol = gRelColumns[i];
        // clean the index columns
        for (uint32_t ci=0; ci<gSchema[i]; ++ci) {
            auto& cCol = cRelCol.columns[ci];
            cCol.transactions.erase(cCol.transactions.begin(),
                    upper_bound(cCol.transactions.begin(), cCol.transactions.end(),
                        f.transactionId,
                        [](const uint64_t target, const ColumnTransaction_t& ct){ return target < ct.first; }
                        ));

        }
        // clean the transactions log
        auto& transLog = gRelations[i].transLog; 
        //cerr << "size bef: " << transLog.size() << endl;
        for (auto it = transLog.begin(), tend=transLog.end(); it!=tend && ((*it)->trans_id <= f.transactionId); ) {
            if ((*it)->aliveTuples == 0 && (*it)->last_del_id <= f.transactionId) { it = transLog.erase(it); tend=transLog.end(); }
            else ++it;
        }

        // delete the transLogTuples
        auto& transLogTuples = gRelations[i].transLogTuples;
        transLogTuples.erase(transLogTuples.begin(), 
                upper_bound(transLogTuples.begin(), transLogTuples.end(), f.transactionId,
                    [](const uint64_t target, const pair<uint64_t, vector<tuple_t>>& o){ return target < o.first; })
                );
        //cerr << "size after: " << transLog.size() << endl;
    }
    /*
    // then delete the transactions from the transaction history 
    auto ub = upper_bound(gTransactionHistory.begin(), 
    gTransactionHistory.end(), 
    f.transactionId,
    [](const uint64_t target, const TransactionStruct& ts){ return target < ts.trans_id; });
    gTransactionHistory.erase(gTransactionHistory.begin(), ub);
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

//static vector<ParseMessageStruct*> gPendingValidationMessages;
static vector<ReceivedMessage*> gPendingValidationMessages;
static std::atomic<uint64_t> gPVMCnt;
/*
   static void processPendingValidationMessagesTask(uint32_t nThreads, uint32_t tid) {
   (void)tid; (void)nThreads;// to avoid unused warning
//cerr << "::: tid " << tid << "new" << endl;
const uint64_t pvmsz = gPendingValidationMessages.size();
for (uint64_t vmi = gPVMCnt++; likely(vmi < pvmsz); vmi=gPVMCnt++) {
//cerr << tid << ":" << vmi << " ";
ReceivedMessage *msg = gPendingValidationMessages[vmi];
processValidationQueries(*reinterpret_cast<const ValidationQueries*>(msg->data.data()), msg->data, tid); 
delete msg;
}
}

static void omp_processPendingValidationMessagesTask(uint32_t nThreads) {
(void)nThreads;// to avoid unused warning
//cerr << "::: tid " << tid << "new" << endl;
const uint64_t pvmsz = gPendingValidationMessages.size();
for (uint64_t vmi = 0; likely(vmi < pvmsz); ++vmi) {
uint32_t tid = omp_get_thread_num();
//cerr << tid << "/" << omp_get_num_threads() << ":" << vmi << endl;
ReceivedMessage *msg = gPendingValidationMessages[vmi];
processValidationQueries(*reinterpret_cast<const ValidationQueries*>(msg->data.data()), msg->data, tid); 
delete msg;
}
}

static void parsePendingValidationMessages(SingleTaskPool &pool, uint32_t nThreads = 4) {
if (unlikely(gPendingValidationMessages.empty())) return;
#ifdef LPDEBUG
auto start = LPTimer.getChrono();
#endif
//cerr << "parsing session: " << gPendingValidationMessages.size() << endl;
(void)nThreads;
pool.startAll(processPendingValidationMessagesTask);
pool.waitAll();

//(void)pool;
//omp_processPendingValidationMessagesTask(nThreads);

gPendingValidationMessages.clear();
gPVMCnt = 0;

#ifdef LPDEBUG
LPTimer.validations += LPTimer.getChrono(start);
#endif
}
 */



void inline initOpenMP(uint32_t nThreads) {
    //omp_set_dynamic(0);           // Explicitly disable dynamic teams
    omp_set_num_threads(nThreads);  // Use 4 threads for all consecutive parallel regions
}

void inline initTBB(uint32_t nThreads) {
   (void)nThreads;
   tbb::task_scheduler_init init(tbb::task_scheduler_init::automatic); 
}


int main(int argc, char**argv) {
    uint64_t numOfThreads = 1;
    if (argc > 1) {
        numOfThreads = strtol(argv[1], NULL, 10);
        if (numOfThreads == LONG_MAX)
            numOfThreads = 1;
    }
    //cerr << "Number of threads: " << numOfThreads << endl;
    Globals.nThreads = numOfThreads;

    initOpenMP(numOfThreads);
    initTBB(numOfThreads);

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
        msgReader = ReaderIOFactory::createAsync(ifs, true);
    } else { 
        msgReader = ReaderIOFactory::createAsync(stdin);
    }

    // do some initial reserves or initializations
    gPendingValidations.reserve(4096); 
    for (uint32_t i=0; i<NUM_RELATIONS; ++i) gTransParseMapPhase[i].reserve(512);
    for (uint32_t i=0; i<NUM_RELATIONS; ++i) {
        gRelations[i].transLog.reserve(512);
        gRelations[i].transLogTuples.reserve(1024);
    }
    // allocate global structures based on thread number
    gStats.reset(new StatStruct[numOfThreads+1]);

    // allocate the workers
    SingleTaskPool workerThreads(numOfThreads, processPendingValidationsTask);
    //SingleTaskPool workerThreads(1, processPendingValidationsTask);
    workerThreads.initThreads();
    // leave two available workes - master - Reader
    MultiTaskPool multiPool(numOfThreads-2);
    //MultiTaskPool multiPool(3);
    multiPool.initThreads();
    multiPool.startAll();

    /*
       tbb::parallel_for((size_t)0, numOfThreads, [=] (size_t i) {
       cerr << "worker " << i << endl;
       });
     */
    try {

        //uint64_t msgs = 0;
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
                        multiPool.addTask(parseValidation, static_cast<void*>(msg)); 
                        //processValidationQueries(*reinterpret_cast<const ValidationQueries*>(msg->data.data()), msg); 

                        break;
                    }
                case MessageHead::Transaction: 
#ifdef LPDEBUG
                    ++gTotalTransactions; 
#endif
                    {Globals.state = GlobalState::TRANSACTION;
                        processTransactionMessage(*reinterpret_cast<const Transaction*>(msg->data.data()), msg); 
                        //delete msg;
                        break;
                    }
                case MessageHead::Flush:  
                    // check if we have pending transactions to be processed
                    multiPool.helpExecution();
                    multiPool.waitAll();
                    //parsePendingValidationMessages(workerThreads, numOfThreads);

                    checkPendingValidations(workerThreads);
                    Globals.state = GlobalState::FLUSH;
                    processFlush(*reinterpret_cast<const Flush*>(msg->data.data()), isTestdriver); 
                    delete msg;
                    break;

                case MessageHead::Forget: 
                    // check if we have pending transactions to be processed
                    multiPool.helpExecution();
                    multiPool.waitAll();
                    //parsePendingValidationMessages(workerThreads, numOfThreads);

                    checkPendingValidations(workerThreads);
                    Globals.state = GlobalState::FORGET;
                    processForget(*reinterpret_cast<const Forget*>(msg->data.data())); 
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
                        cerr << "  :::: " << LPTimer << endl << "total validations: " << gTotalValidations << " trans: " << gTotalTransactions << " tuples: " << gTotalTuples << endl; 
#endif              
                        workerThreads.destroy();
                        multiPool.destroy();
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
            const uint64_t *keys = o.keys;
            for (uint32_t c=0; c<o.rowCount; ++c) *ptr++ = *keys++;
            //gRelTransMutex[o.relationId].lock(); 
            gTransParseMapPhase[o.relationId].emplace_back(t.transactionId, true, o.rowCount, ptr-o.rowCount);
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
            const uint64_t *vptr=o.values;
            uint32_t sz = relCols*o.rowCount;
            for (uint32_t c=0; c<sz; ++c) *tptr++ = *vptr++;
            //gRelTransMutex[o.relationId].lock(); 
            gTransParseMapPhase[o.relationId].emplace_back(t.transactionId, false, o.rowCount, tptr-sz);
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
static std::atomic<uint64_t> gNextIndex;

static void updateRequiredColumns(uint64_t ri, vector<SColType>::iterator colBegin, vector<SColType>::iterator colEnd) {
    // PHASE TWO OF THE ALGORITHM IN THIS STAGE IS TO INCREMENTALLY UPDATE
    // THE INDEXES ONLY FOR THE COLUMNS THAT ARE GOING TO BE REQUESTED IN THE 
    // FOLOWING VALIDATION SESSION - 1st predicates only for now
    //for (SColType& cp : *statCols) cerr << "is Op::Equal " << cp.first << " col: " << cp.second << endl; 
    auto& relation = gRelations[ri];
    auto& relColumns = gRelColumns[ri].columns;

    (void)colBegin; (void)colEnd;
    // for each column to be indexed
#pragma omp parallel for schedule(static, 1) num_threads(2)
    for (uint32_t col=0; col<gSchema[ri]; ++col) {
        //tbb::parallel_for ((uint32_t)0, gSchema[ri], [&] (uint32_t col) {

        //uint32_t rel,col;
        //for (; colBegin!=colEnd; ++colBegin) {
        //    lp::validation::unpackRelCol(colBegin->second, rel, col);
        //cerr << "relation: " << ri << " got rel " << rel << " col " << col << endl;
        auto& colTransactions = relColumns[col].transactions;
        uint64_t updatedUntil = relColumns[col].transTo;

        // Use lower_bound to automatically jump to the transaction to start
        auto transFrom = lower_bound(relation.transLogTuples.begin(), relation.transLogTuples.end(), updatedUntil, TransLogComp);
        // for all the transactions in the relation
        for(auto tEnd=relation.transLogTuples.end(); transFrom!=tEnd; ++transFrom) {
            auto& trp = *transFrom;
            // allocate vectors for the current new transaction to put its data
            colTransactions.emplace_back(trp.first, move(vector<CTransStruct>()));
            colTransactions.back().second.reserve(trp.second.size());
            for (auto tpl : trp.second) {
                colTransactions.back().second.emplace_back(tpl[col], tpl);
            }
            sort(colTransactions.back().second.begin(), colTransactions.back().second.end(), ColTransValueLess);
        }
        if(likely(!relation.transLogTuples.empty()))
            relColumns[col].transTo = max(relation.transLogTuples.back().first+1, updatedUntil);
        //cerr << "col " << col << " ends to " << relColumns[col].transTo << endl;
        //});
    }
}

static void processPendingIndexTask(uint32_t nThreads, uint32_t tid) {
    (void)tid; (void)nThreads;// to avoid unused warning
    //cerr << "::: tid " << tid << "new" << endl;

    for (uint64_t ri = gNextIndex++; likely(ri < NUM_RELATIONS); ri=gNextIndex++) {
        //#pragma omp parallel for schedule(static, 1)  
        //for(uint64_t ri = 0; ri < NUM_RELATIONS; ++ri) {

        //auto colpair = std::equal_range(gStatColumns->begin(), gStatColumns->end(), ri, StatCompRel);
        //auto colBegin = colpair.first, colEnd = colpair.second; 
        auto colBegin = gStatColumns->begin(), colEnd = gStatColumns->end(); 

        // take the vector with the transactions and sort it by transaction id in order to apply them in order
        auto& relTrans = gTransParseMapPhase[ri];
        if (relTrans.empty()) { 
            // TODO - we have to run this regardless of transactions since some
            // columns might have to use previous transactions and be called for the first time
            updateRequiredColumns(ri, colBegin, colEnd);
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
        operations.reserve(4096);
        for (auto& trans : relTrans) {
            if (trans.trans_id != lastTransId) {
                // store the tuples for the last transaction just finished
                if (!operations.empty())
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
                        auto tit = lower_bound(relation.transLog.begin(), relation.transLog.end(), lb->second.first, RTLComp);
                        (*tit)->last_del_id = trans.trans_id;
                        --(*tit)->aliveTuples;

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
        updateRequiredColumns(ri, colBegin, colEnd);
    } // end of while true
}

static inline void checkPendingTransactions(SingleTaskPool& pool) {
#ifdef LPDEBUG
    auto startIndex = LPTimer.getChrono();
#endif

    //cerr << "::: session start ::::" << endl;
    //vector<SColType>* cols = &gStats[0].reqCols;
    /*
       uint64_t totalCols = 0;
    //for (SColType& cp : *cols) cerr << "==: " << cp.first << " col: " << cp.second << endl; 
    for (uint32_t tid=1; tid<Globals.nThreads; ++tid) {
    if (gStats[tid].reqCols.size() > cols->size()) cols = &gStats[tid].reqCols;
    totalCols += gStats[tid].reqCols.size();
    }
    cols->reserve(totalCols);
    for (uint32_t tid=0; tid<Globals.nThreads; ++tid) {
    auto ccols = &gStats[tid].reqCols;
    if (unlikely(ccols == cols || ccols->empty())) continue;
    cols->insert(cols->end(), ccols->begin(), ccols->end());
    ccols->resize(0);
    }

    // add the first column for all relations
    std::sort(cols->begin(), cols->end(), StatComp);
    auto it = std::unique(cols->begin(), cols->end(), StatCompEq);
    cols->erase(it, cols->end());
    //cerr << "unique reqCols: " << cols->size() << endl;
    //for (auto& p : *cols) cerr << " " << p.second;
    //for (SColType& cp : *cols) cerr << "==: " << cp.first << " col: " << cp.second << endl; 
     */
    /*// - insert all the columns for all relations 
      for (uint32_t r=0; r<NUM_RELATIONS; ++r) {
      for (uint32_t c=0; c<gSchema[r]; ++c)
      cols->emplace_back(false, lp::validation::packRelCol(r, c));
      }
     */
    //gStatColumns = cols;
    //for (SColType& cp : *gStatColumns) cerr << "==: " << cp.first << " col: " << cp.second << endl; 
    gNextIndex = 0;
    pool.startAll(processPendingIndexTask);
    pool.waitAll();

    //(void)pool;
    //processPendingIndexTask(4, 0);

    //#pragma omp parallel for schedule(static, 1) num_threads(2)
    for (uint32_t r=0; r<NUM_RELATIONS; ++r) gTransParseMapPhase[r].clear();
    /*
       tbb::parallel_for ((uint32_t)0, NUM_RELATIONS, [&] (uint32_t r) {
       gTransParseMapPhase[r].clear();
       });
     */
    // clear the 1st predicate columns 
    //cols->resize(0);

#ifdef LPDEBUG
    LPTimer.transactionsIndex += LPTimer.getChrono(startIndex);
#endif
}


static uint64_t resIndexOffset = 0;
static std::atomic<uint64_t> gNextPending;

static void checkPendingValidations(SingleTaskPool &pool) {
    if (unlikely(gPendingValidations.empty())) return;

    // check if there is any pending index creation to be made before checking validation
    checkPendingTransactions(pool);

#ifdef LPDEBUG
    auto start = LPTimer.getChrono();
#endif

    resIndexOffset = UINT64_MAX;
    for (auto& pv : gPendingValidations) if (pv.validationId < resIndexOffset) resIndexOffset = pv.validationId;
    auto gPRsz = gPendingResults.size();
    if (gPVunique > gPRsz)
        gPendingResults.resize(gPVunique);
    //memset(gPendingResults.data(), 0, sizeof(PendingResultType)*gPRsz); // TODO - maybe memset better
    std::fill(gPendingResults.begin(), gPendingResults.end(), 0);
    gNextPending.store(0);

    //pool.startAll(processPendingValidationsTask);
    //pool.waitAll();

    (void)pool;
    processPendingValidationsTask(0, 0);

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

typedef CTransStruct TupleType;
typedef vector<TupleType> TupleCont;
//typedef vector<Query::Column>::iterator PredIter;
typedef Query::Column* PredIter;

bool inline isTupleConflict(PredIter cbegin, PredIter cend, TupleType& tup) {
    tuple_t& tuple = tup.tuple;
    //if (v.validationId == 4) cerr << "next tuple: " << tuple << endl;
    bool match=true;
    //for (uint32_t cp=pFrom, sz=q.predicates.size(); cp<sz; ++cp) {
    //auto& c = q.predicates[cp];
    for (auto tbegin = cbegin; tbegin<cend; ++tbegin) {
        //auto& c = q.predicates[cp];
        auto& c = *tbegin;
        // make the actual check
        uint64_t tupleValue = tuple[c.column]; 
        uint64_t queryValue = c.value;
        bool result=false;
        switch (c.op) {
            case Op::Equal: 
                result=(tupleValue==queryValue); 
                break;
            case Op::Less: 
                result=(tupleValue<queryValue); 
                break;
            case Op::LessOrEqual: 
                result=(tupleValue<=queryValue); 
                break;
            case Op::Greater: 
                result=(tupleValue>queryValue); 
                break;
            case Op::GreaterOrEqual: 
                result=(tupleValue>=queryValue); 
                break;
            case Op::NotEqual: 
                result=(tupleValue!=queryValue); 
                break;
        } 
        // there is one predicate not true so this whole query on this relation is false
        if (!result) { return false; }
    } // end of single query predicates
    return match;    
}

bool inline  isTupleRangeConflict(vector<TupleType>::iterator tupFrom, vector<TupleType>::iterator tupTo, PredIter cbegin, PredIter cend) {
    //for(; tupFrom!=tupTo; ++tupFrom) {  
    //    if (isTupleConflict(cbegin, cend, tupFrom)) return true;
    //} // end of all tuples for this transaction
    //return false;
    
    return std::find_if(tupFrom, tupTo, [&](TupleType& tup) { return isTupleConflict(cbegin, cend, tup);}) != tupTo;
}


static bool isTransactionConflict(LPQuery& q, uint64_t tri) {
    auto cbegin = reinterpret_cast<Query::Column*>(q.rawQuery->columns),
         cend = cbegin + q.colCountUniq;
    auto pFirst = *reinterpret_cast<Query::Column*>(q.rawQuery->columns);
        
    auto& relColumns = gRelColumns[q.relationId].columns;
    auto& transactions = relColumns[pFirst.column].transactions;
    
    auto iter = &transactions[tri];
    auto& transValues = iter->second;
    decltype(transValues.begin()) tupFrom, tupTo,
        tBegin = transValues.begin(), tEnd=transValues.end();
    uint32_t pFrom = 1;
    // find the valid tuples using range binary searches based on the first predicate
    switch (pFirst.op) {
        case Op::Equal: 
            {auto tp = std::equal_range(tBegin, tEnd, pFirst.value, ColTransValueLess);
                tupFrom = tp.first; tupTo = tp.second;
                break;}
        case Op::Less: 
            tupFrom = tBegin;                    
            tupTo = std::lower_bound(tBegin, tEnd, pFirst.value, ColTransValueLess);                   
            break;
        case Op::LessOrEqual: 
            tupFrom = tBegin;                    
            tupTo = std::upper_bound(tBegin, tEnd, pFirst.value, ColTransValueLess);                   
            break;
        case Op::Greater: 
            tupFrom = std::upper_bound(tBegin, tEnd, pFirst.value, ColTransValueLess);  
            tupTo = tEnd;                   
            break;
        case Op::GreaterOrEqual: 
            tupFrom = std::lower_bound(tBegin, tEnd, pFirst.value, ColTransValueLess);
            tupTo = tEnd;                   
            break;
        default: 
            tupFrom = tBegin;
            tupTo = tEnd;
            pFrom = 0;
    }

    //cerr << "tup diff " << (tupTo - tupFrom) << endl; 
    //if (std::distance(tupFrom, tupTo) == 0) continue;
    if (std::distance(tupFrom, tupTo) == 0) return false;

    if (pFrom == 1) ++cbegin;
    //if (isTupleRangeConflict(tupFrom, tupTo, cpbegin, cend)) return true;
    if (isTupleRangeConflict(tupFrom, tupTo, cbegin, cend)) return true;
    return false;
}

static bool isValidationConflict(LPValidation& v) {
    // no empty validation has none query - ONLY if we did the satisfiability check before
    //if (unlikely(v.queries.empty())) { cerr << "e "; return false; }

    for (auto& q : v.queries) {
        lp::query::preprocess(q);
        if (!lp::query::satisfiable(q)) continue; // go to the next query

        // protect from the case where there is no single predicate
        //if (q.predicates.empty()) { 
        if (q.colCountUniq == 0) { 
            //cerr << "empty: " << v.validationId << endl; 
            auto& transactionsCheck = gRelations[q.relationId].transLogTuples;
            auto transFromCheck = std::lower_bound(transactionsCheck.begin(), transactionsCheck.end(), v.from, TransLogComp);
            if (transFromCheck == transactionsCheck.end()) { 
                // no transactions exist for this query
                continue;
            } else { 
                // transactions exist for this query so it is a conflict
                return true;
            }; 
        }

        // just find the range of transactions we want in this relation
        auto& transactions = gRelColumns[q.relationId].columns[0].transactions;
        auto transFrom = std::lower_bound(transactions.begin(), transactions.end(), v.from, CTRSLessThan);
        auto transTo = std::upper_bound(transFrom, transactions.end(), v.to, CTRSLessThan);
        //cerr << "after: " << v.from << "-" << v.to << "=" << (transTo-transFrom) << " for col: " << pFirst.column << "-" << pFirst.value << endl;
        
        auto trFidx = std::distance(transactions.begin(), transFrom);
        auto trTidx = std::distance(transactions.begin(), transTo);

        volatile bool conflict = false;
        //for(auto iter=transFrom; iter!=transTo; ++iter) {  
        //for(auto tri=trFidx; trFidx<trTidx; ++trFidx) {  
        tbb::parallel_for ((uint64_t)trFidx, (uint64_t)trTidx, [&] (uint64_t tri) {
            if (isTransactionConflict(q, tri)) conflict = true;
        } // end of all the transactions for this relation for this specific query
        ); // tbb parallel for
        if (conflict) return true;
    }// end for all queries
    return false;
}


static void processPendingValidationsTask(uint32_t nThreads, uint32_t tid) {
    (void)tid; (void)nThreads;// to avoid unused warning

    uint64_t totalPending = gPendingValidations.size();
    uint64_t resPos;
    // get a validation ID - atomic operation
    //for (uint64_t vi = gNextPending++; likely(vi < totalPending); vi=gNextPending++) {
        //#pragma omp parallel for schedule(static, 1)
        //for (uint64_t vi = 0; vi < totalPending; ++vi) {
    tbb::parallel_for ((uint64_t)0, totalPending, [&] (uint64_t vi) {
        auto& v = gPendingValidations[vi];
        resPos = v.validationId - resIndexOffset;
        auto& atoRes = gPendingResults[resPos];
        if(isValidationConflict(v)) { atoRes = true; }
        delete v.rawMsg;
    } // while true take more validations 
    ); // for tbb::parallel_for
}

