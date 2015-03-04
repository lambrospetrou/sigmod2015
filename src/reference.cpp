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


#include <iostream>
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

#include "include/circularfifo_memory_relaxed_aquire_release.hpp"
#include "include/concurrentqueue.h"

#include "include/atomicwrapper.hpp"
#include "include/LPTimer.hpp"
#include "include/BoundedQueue.hpp"
#include "include/BoundedAlloc.hpp"
#include "include/SingleTaskPool.hpp"
#include "include/MultiTaskPool.hpp"
#include "include/MultiTaskPoolConc.hpp"
#include "include/LPSpinLock.hpp"

#include "include/cpp_btree/btree_map.h"

#include "include/ReferenceTypes.hpp"
#include "include/LPQueryTypes.hpp"


//#include "include/tbb/tbb.h"



//---------------------------------------------------------------------------
using namespace std;
using namespace lp;

using namespace memory_relaxed_aquire_release;

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

std::ostream& operator<< (std::ostream& os, const LPQuery& o) {
    os << "{" << o.relationId << "-" << o.columnCount << ":: " << o.predicates << "::" << o.satisfiable << "}";
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

// Custom data structures to hold data
struct CTransStruct {
    uint64_t value;
    tuple_t tuple;
    CTransStruct (uint64_t v, tuple_t t) : value(v), tuple(t) {}
    bool operator< (const CTransStruct& o) { 
        return value < o.value;
    }
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
    std::unique_ptr<uint64_t[]> tuples;
    uint64_t rowCount;
    RelTransLog(uint64_t tid, uint64_t* tpl, uint64_t _rowCount) : trans_id(tid), last_del_id(tid), aliveTuples(_rowCount), rowCount(_rowCount) {
        tuples.reset(tpl);
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
    vector<unique_ptr<RelTransLog>> transLog;
    vector<pair<uint64_t, vector<tuple_t>>> transLogDel;
    unordered_map<uint32_t, pair<uint64_t, uint64_t*>> insertedRows;
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
    bool isDelOp;
    uint32_t rowCount;
    uint64_t *values; // delete op => row keys to delete | insert => tuples

    TRMapPhase(uint64_t tid, bool isdel, uint32_t rows, uint64_t *vals)
        : trans_id(tid), isDelOp(isdel), rowCount(rows), values(vals) {}
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

// TRANSACTION HISTORY STRUCTURES

struct TransOperation {
    uint32_t rel_id;
    std::unique_ptr<uint64_t[]> tuple;
    TransOperation(uint32_t relid, std::unique_ptr<uint64_t[]> t):
        rel_id(relid), tuple(move(t)) {}
};
struct TransactionStruct {
    uint64_t trans_id;
    //vector<TransOperation> operations;
    std::atomic<uint64_t> aliveTuples;
    std::unique_ptr<uint64_t[]> tuples;
    uint64_t rowCount;
    TransactionStruct(uint64_t tid, uint64_t* tpl, uint64_t rowCount) : trans_id(tid) {
        tuples.reset(tpl);
        aliveTuples = rowCount;
    } 
};
struct TSComp_t {
    inline bool ByTrans(const TransactionStruct& l, const TransactionStruct& r) {
        return l.trans_id < r.trans_id;
    }
    inline bool operator() (const std::unique_ptr<TransactionStruct>& l, const std::unique_ptr<TransactionStruct>& r) {
        return l->trans_id < r->trans_id;
    }
    inline bool operator() (const std::unique_ptr<TransactionStruct>& o, uint64_t target) {
        return o->trans_id < target;
    }
} TSComp;
//static vector<std::unique_ptr<TransactionStruct>> gTransactionHistory;
//static std::mutex gTransactionHistoryMutex;
//static LPSpinLock gTransactionHistoryMutex;

static uint32_t NUM_RELATIONS;
static std::unique_ptr<uint32_t[]> gSchema;

///////// AUXILIARY STRUCTURES FOR THE WHOLE PROGRAM

static vector<LPValidation> gPendingValidations;
//typedef atomic_wrapper<bool> PendingResultType;
typedef char PendingResultType;
static vector<PendingResultType> gPendingResults;
static vector<pair<uint64_t,bool>> gQueryResults;
static uint64_t gPVunique;

static std::mutex gPendingValidationsMutex;
//static LPSpinLock gPendingValidationsMutex;


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

struct ReceivedMessage {
    MessageHead head;
    vector<char> data;
};
struct ParseMessageStruct {
    ReceivedMessage *msg;
    uint64_t refId;
    BoundedQueue<ReceivedMessage> *msgQ;
    uint64_t memRefId;
    BoundedAlloc<ParseMessageStruct> *memQ;
};


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

static void processTransactionMessage(const Transaction& t, vector<char>& data);
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
    //cerr << "columns: " << NUM_RELATIONS << endl;
    for(uint32_t ci=0; ci<d.relationCount; ++ci) {
        //gRelColumns[ci].columns.resize(gSchema[ci]);
        gRelColumns[ci].columns.reset(new ColumnStruct[gSchema[ci]]);
        //cerr << " " << gSchema[ci];
    }
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
static void processValidationQueries(const ValidationQueries& v, const vector<char>& vdata, uint64_t tid = 0) {
#ifdef LPDEBUG
    auto start = LPTimer.getChrono();    
#endif
    (void)vdata; (void)tid;


    // TODO - OPTIMIZATION CAN BE DONE IF I JUST COPY THE WHOLE DATA instead of parsing it
    // try to put all the queries into a vector
    vector<LPQuery> queries;
    const char* qreader=v.queries;
    const Query *q;
    for (unsigned int i=0;i<v.queryCount;++i) {
        q=reinterpret_cast<const Query*>(qreader);
        LPQuery nQ;
        //LPQuery nQ(*q);
        //cerr << v.validationId << "====" << v.from << ":" << v.to << nQ << endl;
        //if (!lp::validation::isQueryUnsolvable(nQ)) {
        if (lp::query::parse(q, gSchema[q->relationId], &nQ)) {
            // this is a valid query
            nQ.relationId = q->relationId;
            if (!nQ.predicates.empty()) {
                //std::sort(nQ.predicates.begin(), nQ.predicates.end(), ColumnCompOp);


                // print the proper predicates passed the checks
                //cerr << endl <<  "===== proper predicates" << endl;
                //for (auto& c : nQ.predicates) cerr << c << " " << endl;
                //cerr << "===== new predicates" << endl;
                //lp::query::parse(q, gSchema[q->relationId], &nQ);


                // gather statistics    
                auto& p = nQ.predicates[0];
                uint32_t rc = lp::validation::packRelCol(nQ.relationId, p.column);
                /* 
                   uint32_t rel, col;
                   lp::validation::unpackRelCol(rc, rel, col);
                   if (rel != nQ.relationId || col != p.column) { cerr << "error packing relcol" << endl; exit(1); }
                 */
                //cerr << " " << rc;
                gStats[tid].reqCols.emplace_back((p.op == Op::Equal), rc);
                //cerr << " " << gStats[tid].reqCols.back().second;
            }
            queries.push_back(move(nQ));
        }

        //        queries.push_back(move(LPQuery(*q)));
        qreader+=sizeof(Query)+(sizeof(Query::Column)*q->columnCount);
    }
    //  cerr << v.validationId << "====" << v.from << ":" << v.to << "=" << v.queryCount << "=" << queries << endl;

    gPendingValidationsMutex.lock();
    gPendingValidations.emplace_back(v.validationId, v.from, v.to, move(queries));    
    // update the global pending validations to reflect this new one
    ++gPVunique;
    gPendingValidationsMutex.unlock();

#ifdef LPDEBUG
    LPTimer.validations += LPTimer.getChrono(start);
#endif

    return;
    }

    //---------------------------------------------------------------------------
    static void processFlush(const Flush& f) {
#ifdef LPDEBUG
        auto start = LPTimer.getChrono();
#endif
        static char zero = '0';
        static char one = '1';
        if (!gQueryResults.empty()) {
            uint64_t removed = 0;
            for (auto& vp : gQueryResults) {
                if (vp.first > f.validationId) break;
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

        //cerr << "\nfinished flushing " << f.validationId << endl;
#ifdef LPDEBUG
        LPTimer.flushes += LPTimer.getChrono(start);
#endif
    }
    //---------------------------------------------------------------------------

    static void processForget(const Forget& f) {
#ifdef LPDEBUG
        //cerr << "Forget: " << f.transactionId << endl;
        auto start = LPTimer.getChrono();
#endif
        //(void)f.transactionId;

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

            // delete the transLogDel
            auto& transLogDel = gRelations[i].transLogDel;
            transLogDel.erase(transLogDel.begin(), 
                    upper_bound(transLogDel.begin(), transLogDel.end(), f.transactionId,
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

    typedef CircularFifo<ReceivedMessage*, 5000> LPMsgQ;


    void ReaderTask(LPMsgQ& msgQ) {
#ifdef LPDEBUG
        auto start = LPTimer.getChrono();
#endif
        std::ios_base::sync_with_stdio(false);
        setvbuf ( stdin , NULL , _IOFBF , 1<<22 );
        while (true) {
            // request place from the message queue - it blocks if full
            ReceivedMessage *msg = new ReceivedMessage();
            //auto& head = msg->head;
            //auto& buffer = msg->data;
            /*
            // read the head of the message - type and len
            // Read the message body and cast it to the desired type
            cin.read(reinterpret_cast<char*>(&head),sizeof(head));
            if (!cin) { cerr << "read error" << endl; abort(); } // crude error handling, should never happen

            if (head.type == MessageHead::Done) {
            msgQ.registerEnq(res.refId); 
            // exit the loop since the reader has finished its job
            break;        
            }

            // read the actual message content
            if (head.messageLen > buffer.size()) buffer.resize(head.messageLen);
            cin.read(buffer.data(), head.messageLen);
             */
            auto& head = msg->head;
            auto& msgData = msg->data;
            size_t rd = fread(reinterpret_cast<char*>(&head), sizeof(head), 1, stdin);
            if (rd < 1) { cerr << "read error" << endl; abort(); } // crude error handling, should never happen
#ifdef LPDEBUG // I put the inner timer here to avoid stalls in the msgQ
            auto startInner = LPTimer.getChrono();
#endif

            if (unlikely(head.type == MessageHead::Done)) {
                // exit the loop since the reader has finished its job
                while (!msgQ.push(msg)) { lp_spin_sleep(); }
#ifdef LPDEBUG
                LPTimer.readingTotal += LPTimer.getChrono(start);
#endif
                return;
            }

            // read the actual message content
            msgData.reserve(head.messageLen);
            msgData.resize(head.messageLen);
            rd = fread(reinterpret_cast<char*>(msgData.data()), 1, head.messageLen, stdin);
            // crude error handling, should never happen
            if (rd < head.messageLen) { cerr << "read error" << endl; abort(); }                 
#ifdef LPDEBUG
            LPTimer.reading += LPTimer.getChrono(startInner);
#endif 
            while (!msgQ.push(msg)) { /*cerr << "r" << std::endl;*/ lp_spin_sleep(std::chrono::microseconds(0)); }
        }
        return;
    }

    void inline parseValidation(uint32_t nThreads, uint32_t tid, void *args) {
        (void)tid; (void)nThreads;
        ParseMessageStruct *pvs = static_cast<ParseMessageStruct*>(args);
        processValidationQueries(*reinterpret_cast<const ValidationQueries*>(pvs->msg->data.data()), pvs->msg->data, tid); 
        //pvs->msgQ->registerDeq(pvs->refId);
        //pvs->memQ->free(pvs->memRefId);
        delete pvs->msg;
        delete pvs;
    }


    template<typename T>
        //using ConcQ = tbb::concurrent_queue<T>;
        using ConcQ = moodycamel::ConcurrentQueue<T>;

    static std::atomic<bool> gPendingValQueriesFinished;
    static ConcQ<ReceivedMessage*> gPendingValQueries;
    static std::atomic<uint64_t> gPendingValQueriesCount;

    static inline void processValidationMessages(uint32_t nThreads, uint32_t tid) {
        (void)nThreads; (void)tid;
        //cerr << "tid: " << tid << endl;
        vector<ReceivedMessage*> msgs; msgs.resize(100);
        //for (; !gPendingValQueriesFinished || (gPendingValQueriesFinished && !gPendingValQueries.empty()); ) {
        while(!gPendingValQueriesFinished 
                || (gPendingValQueriesFinished 
                    && gPendingValQueriesCount.load(std::memory_order_acquire) > 0)) {

            //ReceivedMessage *msg;
            //while (!gPendingValQueries.try_dequeue(msg)) {
            std::size_t res;
            while ((res = gPendingValQueries.try_dequeue_bulk(msgs.data(), 100)) != 0) {
                gPendingValQueriesCount.fetch_sub(res, std::memory_order_release);
                for (std::size_t i=0; i<res; ++i) {
                    auto msg = msgs[i];
                    processValidationQueries(*reinterpret_cast<const ValidationQueries*>(msg->data.data()), msg->data, tid); 
                    delete msg;
                }
            }
            // failed to get jobs so check again
            if (gPendingValQueriesFinished) return; // no more messages are coming
            cerr << "w ";
            lp_spin_sleep(std::chrono::microseconds((tid&1)?25:100));
            //lp_spin_sleep(std::chrono::microseconds(0));
        }// end of while there are messages
        }

        void inline parseTransactionPH1(uint32_t nThreads, uint32_t tid, void *args) {
            (void)tid; (void)nThreads;
            ParseMessageStruct *pvs = static_cast<ParseMessageStruct*>(args);
            processTransactionMessage(*reinterpret_cast<const Transaction*>(pvs->msg->data.data()), pvs->msg->data); 
            //pvs->msgQ->registerDeq(pvs->refId);
            //pvs->memQ->free(pvs->memRefId);
            delete pvs->msg;
            delete pvs;
        }


#ifdef LPDEBUG
        static uint64_t gTotalTransactions = 0, gTotalTuples = 0, gTotalValidations = 0;
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

            //const uint64_t MessageQSize = 500;
            //BoundedQueue<ReceivedMessage> msgQ(MessageQSize);
            //BoundedAlloc<ParseMessageStruct> memQ(MessageQSize);

            LPMsgQ msgQ;

            std::thread readerTask(ReaderTask, std::ref(msgQ));

            gPendingValQueriesFinished = false;

            //SingleTaskPool workerThreads(numOfThreads, processValidationMessages);
            SingleTaskPool workerThreads(numOfThreads, processPendingValidationsTask);
            workerThreads.initThreads();
            //workerThreads.startAll();
            // leave two available workes - master - msgQ
            MultiTaskPool multiPool(numOfThreads-2);
            //MultiTaskPool multiPool(1);
            multiPool.initThreads();
            multiPool.startAll();

            // allocate global structures based on thread number
            gStats.reset(new StatStruct[numOfThreads+1]);

            vector<ReceivedMessage*> msgs;

            try {

                //uint64_t msgs = 0;
                while (true) {
                    /*
#ifdef LPDEBUG
auto start = LPTimer.getChrono();
#endif
ReceivedMessage *msg = new ReceivedMessage();
auto& head = msg->head;
auto& msgData = msg->data;
                    // read the head of the message - type and len
                    // Read the message body and cast it to the desired type
                    //cin.read(reinterpret_cast<char*>(&head),sizeof(head));
                    //if (!cin) { cerr << "read error" << endl; abort(); } // crude error handling, should never happen
                    size_t rd = fread(reinterpret_cast<char*>(&head), sizeof(head), 1, stdin);
                    if (rd < 1) { cerr << "read error" << endl; abort(); } // crude error handling, should never happen
                    // read the actual message content
                    msgData.reserve(head.messageLen);
                    msgData.resize(head.messageLen);
                    rd = fread(reinterpret_cast<char*>(msgData.data()), 1, head.messageLen, stdin);
                    if (rd < head.messageLen) { cerr << "read error" << endl; abort(); } // crude error handling, should never happen
#ifdef LPDEBUG
LPTimer.readingTotal += LPTimer.getChrono(start);
#endif 
                     */

                    ReceivedMessage *msg;
                    while (!msgQ.pop(msg)) {/*cerr<<"m ";*/ lp_spin_sleep(std::chrono::microseconds(50));}
                    auto& head = msg->head;
                    //auto& msgData = msg->data;
                    // Retrieve the message
                    //cerr << "try for incoming" << endl;
                    //auto res = msgQ.reqNextDeq();
                    //cerr << "deq id: " << res.refId << endl;
                    //ReceivedMessage& msg = *res.value;
                    //auto& head = msg.head;
                    //cerr << "incoming: " << head.type << " : " << head.messageLen << " data: " << msgData.size() << endl;
                    // And interpret it
                    switch (head.type) {
                        case MessageHead::ValidationQueries: 
                            {    Globals.state = GlobalState::VALIDATION;
                                //processValidationQueries(*reinterpret_cast<const ValidationQueries*>(msg.data.data()), msg.data); 
                                //msgQ.registerDeq(res.refId);
#ifdef LPDEBUG
                                ++gTotalValidations; // this is just to count the total validations....not really needed!
#endif
                                //++gPendingValQueriesCount;
                                /*
                                   msgs.push_back(msg);
                                   if (msgs.size() > 50) {
                                   gPendingValQueries.enqueue_bulk(msgs.data(), msgs.size());
                                   gPendingValQueriesCount.fetch_add(msgs.size(), std::memory_order_release);
                                   msgs.resize(0);
                                   }
                                 */
                                //cerr << "inserted msg to conc queue" << endl;
                                ParseMessageStruct *pvs = new ParseMessageStruct();
                                pvs->msg = msg;
                                multiPool.addTask(parseValidation, static_cast<void*>(pvs)); 

                                break;
                            }
                        case MessageHead::Transaction: 
#ifdef LPDEBUG
                            ++gTotalTransactions; 
#endif
                            {Globals.state = GlobalState::TRANSACTION;
                                //processTransactionMessage(*reinterpret_cast<const Transaction*>(msg.data.data()), msg.data); 
                                processTransactionMessage(*reinterpret_cast<const Transaction*>(msg->data.data()), msg->data); 
                                delete msg;
                                //ParseMessageStruct *pvs = new ParseMessageStruct();
                                //pvs->msg = msg;
                                //parseTransactionPH1(numOfThreads, numOfThreads, pvs);
                                //msgQ.registerDeq(res.refId);
                                /* 
                                   BoundedAlloc<ParseMessageStruct>::BAResult& mem = memQ.malloc();
                                   ParseMessageStruct *pvs = mem.value;
                                   pvs->msgQ = &msgQ;
                                   pvs->refId = res.refId;
                                   pvs->memRefId = mem.refId;
                                   pvs->memQ = &memQ;
                                   pvs->msg = &msg;
                                //parseTransactionPH1(numOfThreads, numOfThreads, pvs); 
                                multiPool.addTask(parseTransactionPH1, static_cast<void*>(pvs)); 
                                 */
                                break;
                            }
                        case MessageHead::Flush:  
                            // check if we have pending transactions to be processed
                            multiPool.helpExecution();
                            /*
                               if (!msgs.empty()) {
                               gPendingValQueries.enqueue_bulk(msgs.data(), msgs.size());
                               gPendingValQueriesCount.fetch_add(msgs.size(), std::memory_order_release);
                               msgs.resize(0);
                               }

                               gPendingValQueriesFinished = true;
                               workerThreads.waitAll();
                             */

                            //cerr << "left msgs: " << gPendingValQueriesCount << endl;

                            multiPool.waitAll();
                            checkPendingValidations(workerThreads);
                            Globals.state = GlobalState::FLUSH;
                            processFlush(*reinterpret_cast<const Flush*>(msg->data.data())); 
                            //msgQ.registerDeq(res.refId);
                            delete msg;
                            /*
                               gPendingValQueriesFinished = false; 
                               workerThreads.startAll(processValidationMessages);
                             */
                            break;

                        case MessageHead::Forget: 
                            // check if we have pending transactions to be processed
                            /*
                               if (!msgs.empty()) {
                               gPendingValQueries.enqueue_bulk(msgs.data(), msgs.size());
                               gPendingValQueriesCount.fetch_add(msgs.size(), std::memory_order_release);
                               msgs.resize(0);
                               }
                               gPendingValQueriesFinished = true;
                               workerThreads.waitAll();
                             */
                            //cerr << "left msgs: " << gPendingValQueriesCount << endl;

                            multiPool.helpExecution();
                            multiPool.waitAll();
                            checkPendingValidations(workerThreads);
                            Globals.state = GlobalState::FORGET;
                            processForget(*reinterpret_cast<const Forget*>(msg->data.data())); 
                            //msgQ.registerDeq(res.refId);
                            delete msg;
                            /*
                               gPendingValQueriesFinished = false; 
                               workerThreads.startAll(processValidationMessages);
                             */
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
                                //gPendingValQueriesFinished = true;
                                //workerThreads.waitAll();
                                //cerr << "left msgs: " << gPendingValQueriesCount << endl;

                                readerTask.join();
                                workerThreads.destroy();
                                multiPool.destroy();
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

        static void processTransactionMessage(const Transaction& t, vector<char>& data) {
#ifdef LPDEBUG
            auto start = LPTimer.getChrono();
#endif
            (void)data;

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
                    gRelTransMutex[o.relationId].lock(); 
                    gTransParseMapPhase[o.relationId].emplace_back(t.transactionId, true, o.rowCount, ptr-o.rowCount);
                    gRelTransMutex[o.relationId].unlock(); 
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
                    gRelTransMutex[o.relationId].lock(); 
                    gTransParseMapPhase[o.relationId].emplace_back(t.transactionId, false, o.rowCount, tptr-sz);
                    gRelTransMutex[o.relationId].unlock(); 
#ifdef LPDEBUG
                    ++gTotalTuples;
#endif
                }// end of lock_guard
                // advance to next Relation insertions
                reader+=sizeof(TransactionOperationInsert)+(sizeof(uint64_t)*o.rowCount*relCols);
            }
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

            // for each column to be indexed
            uint32_t rel,col;
            for (; colBegin!=colEnd; ++colBegin) {
                lp::validation::unpackRelCol(colBegin->second, rel, col);
                //cerr << "relation: " << ri << " got rel " << rel << " col " << col << endl;
                auto& colTransactions = relColumns[col].transactions;
                uint64_t updatedUntil = relColumns[col].transTo;

                // Use lower_bound to automatically jump to the transaction to start
                auto transFrom = lower_bound(relation.transLogDel.begin(), relation.transLogDel.end(), updatedUntil, TransLogComp);
                // for all the transactions in the relation
                for(auto tEnd=relation.transLogDel.end(); transFrom!=tEnd; ++transFrom) {
                    auto& trp = *transFrom;
                    // allocate vectors for the current new transaction to put its data
                    colTransactions.emplace_back(trp.first, move(vector<CTransStruct>()));
                    colTransactions.back().second.reserve(trp.second.size());
                    for (auto tpl : trp.second) {
                        colTransactions.back().second.emplace_back(tpl[col], tpl);
                    }
                    sort(colTransactions.back().second.begin(), colTransactions.back().second.end(), ColTransValueLess);
                }
                if(!relation.transLogDel.empty())
                    relColumns[col].transTo = max(relation.transLogDel.back().first+1, updatedUntil);
                //cerr << "col " << col << " ends to " << relColumns[col].transTo << endl;
            }
        }

        static void processPendingIndexTask(uint32_t nThreads, uint32_t tid) {
            (void)tid; (void)nThreads;// to avoid unused warning
            //cerr << "::: tid " << tid << "new" << endl;

            for (;true;) {
                uint64_t ri = gNextIndex++;
                if (ri >= NUM_RELATIONS) return;

                auto colpair = std::equal_range(gStatColumns->begin(), gStatColumns->end(), ri, StatCompRel);
                auto colBegin = colpair.first, colEnd = colpair.second; 

                // take the vector with the transactions and sort it by transaction id in order to apply them in order
                auto& relTrans = gTransParseMapPhase[ri];
                if (relTrans.empty()) { 
                    // TODO - we have to run this regardless of transactions since some
                    // columns might have to use previous transactions and be called for the first time
                    updateRequiredColumns(ri, colBegin, colEnd);
                    continue; 
                }

                //cerr << "tid " << tid << " got " << ri << " = " << relTrans.size() << endl;

                std::sort(relTrans.begin(), relTrans.end(), TRMapPhaseByTrans);

                auto& relation = gRelations[ri];
                //auto& relColumns = gRelColumns[ri].columns;
                uint32_t relCols = gSchema[ri];

                uint64_t lastTransId = relTrans[0].trans_id;

                // for each transaction regarding this relation
                vector<tuple_t> operations;
                for (auto& trans : relTrans) {
                    if (trans.trans_id != lastTransId) {
                        // store the tuples for the last transaction just finished
                        if (!operations.empty())
                            relation.transLogDel.emplace_back(lastTransId, operations);
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
                                tuple_t tpl = lb->second.second;
                                operations.push_back(tpl);

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
                        // TODO - THIS HAS TO BE IN ORDER - each transaction will have its own transaction history from now on
                        relation.transLog.emplace_back(new RelTransLog(trans.trans_id, trans.values, trans.rowCount));
                    }
                }
                // store the last transaction data
                // store the operations for the last transaction
                if (!operations.empty())
                    relation.transLogDel.emplace_back(lastTransId, move(operations));

                // update with new transactions
                updateRequiredColumns(ri, colBegin, colEnd);
            } // end of while true
        }

        static inline void checkPendingTransactions(SingleTaskPool& pool) {
#ifdef LPDEBUG
            auto startIndex = LPTimer.getChrono();
#endif

            //cerr << "::: session start ::::" << endl;
            vector<SColType>* cols = &gStats[0].reqCols;
            uint64_t totalCols = 0;
            //for (SColType& cp : *cols) cerr << "==: " << cp.first << " col: " << cp.second << endl; 
            for (uint32_t tid=1; tid<Globals.nThreads; tid++) {
                if (gStats[tid].reqCols.size() > cols->size()) cols = &gStats[tid].reqCols;
                totalCols += gStats[tid].reqCols.size();
            }
            cols->reserve(totalCols);
            for (uint32_t tid=0; tid<Globals.nThreads; tid++) {
                auto ccols = &gStats[tid].reqCols;
                if (ccols == cols) continue;
                //copy(ccols->begin(), ccols->end(), back_inserter(*cols));
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

            gStatColumns = cols;
            //for (SColType& cp : *gStatColumns) cerr << "==: " << cp.first << " col: " << cp.second << endl; 

            gNextIndex = 0;
            pool.startAll(processPendingIndexTask);
            pool.waitAll();
            for (uint32_t r=0; r<NUM_RELATIONS; ++r) gTransParseMapPhase[r].clear();

            // clear the 1st predicate columns 
            cols->resize(0);

#ifdef LPDEBUG
            LPTimer.transactionsIndex += LPTimer.getChrono(startIndex);
#endif
        }


        static uint64_t resIndexOffset = 0;
        static std::atomic<uint64_t> gNextPending;

        static void checkPendingValidations(SingleTaskPool &pool) {
            if (gPendingValidations.empty()) return;

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

            pool.startAll(processPendingValidationsTask);
            pool.waitAll();

            // update the results - you can get the validation id by adding resIndexOffset to the position
            for (uint64_t i=0, valId=resIndexOffset; i<gPVunique; ++i, ++valId) { 
                //gQueryResults.push_back(move(make_pair(valId, gPendingResults[i])));
                gQueryResults.emplace_back(valId, gPendingResults[i]);
            }
            gPendingValidations.clear();
            gPVunique = 0;
#ifdef LPDEBUG
            LPTimer.validationsProcessing += LPTimer.getChrono(start);
#endif
        }


        static void processPendingValidationsTask(uint32_t nThreads, uint32_t tid) {
            (void)tid; (void)nThreads;// to avoid unused warning

            uint64_t totalPending = gPendingValidations.size();
            uint64_t resPos, vi;
            for (;true;) {

                // get a validation ID - atomic operation
                vi = gNextPending++;
                if (vi >= totalPending) { /*cerr << "exiting" << endl;*/  return; } // all pending finished

                auto& v = gPendingValidations[vi];
                resPos = v.validationId - resIndexOffset;
                auto& atoRes = gPendingResults[resPos];
                // check if someone else found a conflict already for this validation ID
                if (atoRes) continue;

                /*
                   for (auto it = v.queries.begin(); it!=v.queries.end();) {
                   std::sort(it->predicates.begin(), it->predicates.end(), );

                   it->predicates.resize(std::distance(it->predicates.begin(), std::unique(it->predicates.begin(), it->predicates.end())));
                   if (lp::validation::isQueryUnsolvable(*it)) {
                   it = v.queries.erase(it);
                   } else {
                   ++it;
                   }
                   }
                 */
                if (v.queries.empty()) { continue; }
                // sort the queries based on the number of the columns needed to check
                // small queries first in order to try finding a solution faster
                std::sort(v.queries.begin(), v.queries.end(), LPQueryCompSizeLess);

                bool conflict = false, otherFinishedThis = false;
                // for each query in this validation         
                for (auto& q : v.queries) {
                    if (atoRes) { otherFinishedThis = true; /*cerr << "h" << endl;*/ break; }

                    // protect from the case where there is no single predicate
                    if (q.predicates.empty()) { 
                        //cerr << "empty: " << v.validationId << endl; 
                        auto& transactionsCheck = gRelations[q.relationId].transLogDel;
                        auto transFromCheck = std::lower_bound(transactionsCheck.begin(), transactionsCheck.end(), v.from, TransLogComp);
                        if (transFromCheck == transactionsCheck.end()) { 
                            // no transactions exist for this query
                            continue;
                        } else { 
                            // transactions exist for this query
                            conflict=true;  break; 
                        }; 
                    }

                    auto& relColumns = gRelColumns[q.relationId].columns;

                    // IMPORTANT!!! - sort them in order to have the equality checks first -  done earlier now
                    //std::sort(q.predicates.begin(), q.predicates.end(), LPQuery::QCSortOp);

                    auto& pFirst = q.predicates[0];
                    auto& transactions = relColumns[pFirst.column].transactions;
                    auto transFrom = std::lower_bound(transactions.begin(), transactions.end(), v.from, CTRSLessThan);
                    auto transTo = std::upper_bound(transFrom, transactions.end(), v.to, CTRSLessThan);

                    //cerr << "after: " << v.from << "-" << v.to << "=" << (transTo-transFrom) << " for col: " << pFirst.column << "-" << pFirst.value << endl;

                    uint32_t pFrom = 0;
                    for(auto iter=transFrom; iter!=transTo; ++iter) {  
                        auto& transValues = iter->second;
                        decltype(transValues.begin()) tupFrom, tupTo, 
                            tBegin = transValues.begin(), tEnd=transValues.end();
                        // find the valid tuples using range binary searches based on the first predicate
                        switch (pFirst.op) {
                            case Op::Equal: 
                                {auto tp = std::equal_range(tBegin, tEnd, pFirst.value, ColTransValueLess);
                                    tupFrom = tp.first; tupTo = tp.second;
                                    pFrom = 1; break;}
                            case Op::Less: 
                                tupFrom = tBegin;                    
                                tupTo = std::lower_bound(tBegin, tEnd, pFirst.value, ColTransValueLess);                   
                                pFrom = 1; break;
                            case Op::LessOrEqual: 
                                tupFrom = tBegin;                    
                                tupTo = std::upper_bound(tBegin, tEnd, pFirst.value, ColTransValueLess);                   
                                pFrom = 1; break;
                            case Op::Greater: 
                                tupFrom = std::upper_bound(tBegin, tEnd, pFirst.value, ColTransValueLess);  
                                tupTo = tEnd;                   
                                pFrom = 1; break;
                            case Op::GreaterOrEqual: 
                                tupFrom = std::lower_bound(tBegin, tEnd, pFirst.value, ColTransValueLess);
                                tupTo = tEnd;                   
                                pFrom = 1; break;
                            default: 
                                tupFrom = tBegin;
                                tupTo = tEnd;
                                pFrom = 0;
                        }

                        //cerr << "tup diff " << (tupTo - tupFrom) << endl; 
                        if (tupTo == tupFrom) continue;

                        for(; tupFrom!=tupTo; ++tupFrom) {  
                            tuple_t& tuple = tupFrom->tuple;
                            //if (v.validationId == 4) cerr << "next tuple: " << tuple << endl;
                            bool match=true;
                            for (uint32_t cp=pFrom, sz=q.predicates.size(); cp<sz; ++cp) {
                                auto& c = q.predicates[cp];
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
                                if (!result) { match=false; break; }
                            } // end of single query predicates
                            //cerr << "match: " << match << " conflict: " << conflict << endl;
                            if (match) { conflict=true; goto CONFLICT;  break; }    
                        } // end of all tuples for this transaction
                        //if (conflict) { break; }
                    } // end of all the transactions for this relation for this specific query
                    //if (conflict) break;
                }// end for all queries
CONFLICT:
                // update the pending results to help other threads skip this validation 
                // if it has other parts
                if (conflict && !otherFinishedThis)  { /*cerr<< "c: " << v.validationId << endl;*/  atoRes = true;}
            } // while true take more validations 
        }

