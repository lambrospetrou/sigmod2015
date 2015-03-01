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
#include <iostream>
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

#include "include/ReferenceTypes.hpp"
#include "include/LPQueryTypes.hpp"

#include "include/atomicwrapper.hpp"
#include "include/LPTimer.hpp"
#include "include/BoundedQueue.hpp"
#include "include/BoundedAlloc.hpp"
#include "include/SingleTaskPool.hpp"
#include "include/MultiTaskPool.hpp"
#include "include/LPSpinLock.hpp"

#include "include/cpp_btree/btree_map.h"

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

std::ostream& operator<< (std::ostream& os, const LPQuery& o) {
    os << "{" << o.relationId << "-" << o.columnCount << ":: " << o.predicates << "::" << o.satisfiable << "}";
    return os;
}
std::ostream& operator<< (std::ostream& os, const Query::Column& o) {
    os << "[" << o.column << ":" << o.op << ":" << o.value << "]";
    return os;
}


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

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

typedef uint64_t* tuple_t;

// Custom data structures to hold data
struct CTransStruct {
    uint64_t value;
    tuple_t tuple;
    CTransStruct (uint64_t v, tuple_t t) : value(v), tuple(t) {}
    bool operator< (const CTransStruct& o) { 
        return value < o.value;
    }
    static bool CompValOnly(const CTransStruct& l, const CTransStruct& r) {
        return l.value < r.value;
    }
};
std::ostream& operator<< (std::ostream& os, const CTransStruct& o) {
    os << "{" << "-" << o.value << "}";
    return os;
}
struct CTRSValueLessThan_t {
    bool operator() (const CTransStruct& o, uint64_t target) {
        return o.value < target;
    }
    bool operator() (uint64_t target, const CTransStruct& o) {
        return target < o.value;
    }
} CTRSValueLessThan;

typedef pair<uint64_t, vector<CTransStruct>> ColumnTransaction_t;
struct ColumnStruct {
    // the trans_id the transactions are updated to inclusive
    uint32_t transTo;
    vector<ColumnTransaction_t> transactions;
};
struct CTRSLessThan_t {
    bool operator() (const ColumnTransaction_t& left, const ColumnTransaction_t& right) {
        return left.first < right.first;
    }
    bool operator() (const ColumnTransaction_t& o, uint64_t target) {
        return o.first < target;
    }
    bool operator() (uint64_t target, const ColumnTransaction_t& o) {
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
    static bool ByTrans(const RelTransLog& l, const RelTransLog& r) {
        return l.trans_id < r.trans_id;
    }
};
struct RTLComp_t {
    bool operator() (const RelTransLog& l, const RelTransLog& r) {
        return l.trans_id < r.trans_id;
    }
    bool operator() (const RelTransLog& o, uint64_t target) {
        return o.trans_id < target;
    }
    bool operator() (const unique_ptr<RelTransLog>& l, const unique_ptr<RelTransLog>& r) {
        return l->trans_id < r->trans_id;
    }
    bool operator() (const unique_ptr<RelTransLog>& o, uint64_t target) {
        return o->trans_id < target;
    }
} RTLComp;
// The general structure for each relation
struct RelationStruct {
    vector<unique_ptr<RelTransLog>> transLog;
    vector<pair<uint64_t, vector<tuple_t>>> transLogDel;
    unordered_map<uint32_t, pair<uint64_t, uint64_t*>> insertedRows;
};

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

    static bool ByTrans(const TransactionStruct& l, const TransactionStruct& r) {
        return l.trans_id < r.trans_id;
    }
};
struct TSComp_t {
    bool operator() (const std::unique_ptr<TransactionStruct>& l, const std::unique_ptr<TransactionStruct>& r) {
        return l->trans_id < r->trans_id;
    }
    bool operator() (const std::unique_ptr<TransactionStruct>& o, uint64_t target) {
        return o->trans_id < target;
    }
} TSComp;
static vector<std::unique_ptr<TransactionStruct>> gTransactionHistory;
static std::mutex gTransactionHistoryMutex;
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
typedef pair<bool, uint32_t> SColType;
struct StatStruct {
    // columns info that appear as 1st predicates - bool=True means equality, False anything else
    vector<SColType> reqCols; 
};
struct StatComp_t {
    bool operator() (const SColType& l, const SColType& r) {
        return l.second < r.second;
    }
    bool operator() (const SColType& l, uint32_t target) {
        return l.second < target;
    }
    /*
    bool operator() (const SColType& l, const SColType& r) {
        if (l.first && !r.first) return true;
        else if (r.first && !l.first) return false;
        return l.second < r.second;
    }
    */
} StatComp;
static unique_ptr<StatStruct[]> gStats;
static vector<SColType> *gStatColumns;

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
        LPQuery nQ(*q);
        //cerr << v.validationId << "====" << v.from << ":" << v.to << nQ << endl;
        if (!lp::validation::isQueryUnsolvable(nQ)) {
            // this is a valid query
            if (!nQ.predicates.empty()) {
                // gather statistics    
                auto& p = nQ.predicates[0];
                //uint32_t rc = ((nQ.relationId & MSK_L16) << 16) | (p.column & MSK_L16);
                uint32_t rc = lp::validation::packRelCol(nQ.relationId, p.column);
                //cerr << " " << rc;
                gStats[tid].reqCols.push_back(move(make_pair((p.op == lp::LPOps::Equal), rc)));
                //cerr << " " << gStats[tid].reqCols.back().second;
            }
            queries.push_back(move(nQ));
        }

//        queries.push_back(move(LPQuery(*q)));
        qreader+=sizeof(Query)+(sizeof(Query::Column)*q->columnCount);
    }
    //  cerr << v.validationId << "====" << v.from << ":" << v.to << "=" << v.queryCount << "=" << queries << endl;
    {
        //std::lock_guard<std::mutex> lk(gPendingValidationsMutex);
        gPendingValidationsMutex.lock();
        gPendingValidations.emplace_back(v.validationId, v.from, v.to, move(queries));    
        // update the global pending validations to reflect this new one
        ++gPVunique;
        gPendingValidationsMutex.unlock();;
    }
#ifdef LPDEBUG
    LPTimer.validations += LPTimer.getChrono(start);
#endif
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
void ReaderTask(BoundedQueue<ReceivedMessage>& msgQ) {
#ifdef LPDEBUG
    auto start = LPTimer.getChrono();
#endif
    while (true) {
        // request place from the message queue - it blocks if full
        auto res = msgQ.reqNextEnq();
        ReceivedMessage& msg = *res.value;
        auto& head = msg.head;
        auto& buffer = msg.data;
        // read the head of the message - type and len
        // Read the message body and cast it to the desired type
        cin.read(reinterpret_cast<char*>(&head),sizeof(head));
#ifdef LPDEBUG // I put the inner timer here to avoid stalls in the msgQ
        auto startInner = LPTimer.getChrono();
#endif
        if (!cin) { cerr << "read error" << endl; abort(); } // crude error handling, should never happen

        if (head.type == MessageHead::Done) {
            msgQ.registerEnq(res.refId); 
            // exit the loop since the reader has finished its job
            break;        
        }

        // read the actual message content
        if (head.messageLen > buffer.size()) buffer.resize(head.messageLen);
        cin.read(buffer.data(), head.messageLen);

#ifdef LPDEBUG
        LPTimer.reading += LPTimer.getChrono(startInner);
#endif 
        msgQ.registerEnq(res.refId);
    }
#ifdef LPDEBUG
    LPTimer.readingTotal += LPTimer.getChrono(start);
#endif 
    return;
}

void inline parseValidation(uint32_t nThreads, uint32_t tid, void *args) {
    (void)tid; (void)nThreads;
    ParseMessageStruct *pvs = static_cast<ParseMessageStruct*>(args);
    processValidationQueries(*reinterpret_cast<const ValidationQueries*>(pvs->msg->data.data()), pvs->msg->data, tid); 
    pvs->msgQ->registerDeq(pvs->refId);
    pvs->memQ->free(pvs->memRefId);
}


void inline parseTransactionPH1(uint32_t nThreads, uint32_t tid, void *args) {
    (void)tid; (void)nThreads;
    ParseMessageStruct *pvs = static_cast<ParseMessageStruct*>(args);
    processTransactionMessage(*reinterpret_cast<const Transaction*>(pvs->msg->data.data()), pvs->msg->data); 
    pvs->msgQ->registerDeq(pvs->refId);
    pvs->memQ->free(pvs->memRefId);
}


#ifdef LPDEBUG
static uint64_t gTotalTransactions = 0, gTotalTuples = 0;
#endif

int main(int argc, char**argv) {
    uint64_t numOfThreads = 1;
    if (argc > 1) {
        numOfThreads = strtol(argv[1], NULL, 10);
        if (numOfThreads == LONG_MAX)
            numOfThreads = 1;
    }
    cerr << "Number of threads: " << numOfThreads << endl;
    Globals.nThreads = numOfThreads;

    uint64_t MessageQSize = 500;
    BoundedQueue<ReceivedMessage> msgQ(MessageQSize);
    BoundedAlloc<ParseMessageStruct> memQ(MessageQSize);

    std::thread readerTask(ReaderTask, std::ref(msgQ));

    //SingleTaskPool workerThreads(numOfThreads, processPendingValidationsTask);
    SingleTaskPool workerThreads(1, processPendingValidationsTask);
    workerThreads.initThreads();

    // leave two available workes - master - msgQ
    MultiTaskPool multiPool(1);
    multiPool.initThreads();
    multiPool.startAll();

    // allocate global structures based on thread number
    gStats.reset(new StatStruct[numOfThreads+1]);

    uint64_t gTotalValidations = 0;

    try {
        //uint64_t msgs = 0;
        while (true) {

            // Retrieve the message
            //cerr << "try for incoming" << endl;
            auto res = msgQ.reqNextDeq();
            //cerr << "deq id: " << res.refId << endl;
            ReceivedMessage& msg = *res.value;
            auto& head = msg.head;
            //cerr << "incoming: " << head.type << " =" << msgs++ << endl;
            // And interpret it
            switch (head.type) {
                case MessageHead::ValidationQueries: 
                    {    Globals.state = GlobalState::VALIDATION;
                        //processValidationQueries(*reinterpret_cast<const ValidationQueries*>(msg.data.data()), msg.data); 
#ifdef LPDEBUG
                        ++gTotalValidations; // this is just to count the total validations....not really needed!
#endif
                        //ParseValidationStruct *pvs = new ParseValidationStruct();
                        BoundedAlloc<ParseMessageStruct>::BAResult& mem = memQ.malloc();
                        ParseMessageStruct *pvs = mem.value;
                        pvs->msgQ = &msgQ;
                        pvs->refId = res.refId;
                        pvs->memRefId = mem.refId;
                        pvs->memQ = &memQ;
                        pvs->msg = &msg;
                        multiPool.addTask(parseValidation, static_cast<void*>(pvs)); 
                        break;
                    }
                case MessageHead::Transaction: 
#ifdef LPDEBUG
                    ++gTotalTransactions; 
#endif
                    {Globals.state = GlobalState::TRANSACTION;
                    processTransactionMessage(*reinterpret_cast<const Transaction*>(msg.data.data()), msg.data); 
                    msgQ.registerDeq(res.refId);
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
                    //multiPool.helpExecution();
                    multiPool.waitAll();
                    checkPendingValidations(workerThreads);
                    Globals.state = GlobalState::FLUSH;
                    processFlush(*reinterpret_cast<const Flush*>(msg.data.data())); 
                    msgQ.registerDeq(res.refId);
                    break;

                case MessageHead::Forget: 
                    // check if we have pending transactions to be processed
                    //multiPool.helpExecution();
                    multiPool.waitAll();
                    checkPendingValidations(workerThreads);
                    Globals.state = GlobalState::FORGET;
                    processForget(*reinterpret_cast<const Forget*>(msg.data.data())); 
                    msgQ.registerDeq(res.refId);
                    break;

                case MessageHead::DefineSchema: 
                    Globals.state = GlobalState::SCHEMA;
                    processDefineSchema(*reinterpret_cast<const DefineSchema*>(msg.data.data()));
                    msgQ.registerDeq(res.refId);
                    break;
                case MessageHead::Done: 
                    {
#ifdef LPDEBUG
                        cerr << "  :::: " << LPTimer << endl << "total validations: " << gTotalValidations << " trans: " << gTotalTransactions << " tuples: " << gTotalTuples << endl; 
#endif              
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

bool TRMapPhaseByTrans(const TRMapPhase& l, const TRMapPhase& r) {
    //return l.trans_id < r.trans_id;
    if (l.trans_id < r.trans_id) return true;
    else if (r.trans_id < l.trans_id) return false;
    else return l.isDelOp && !r.isDelOp;
}

static void processPendingIndexTask(uint32_t nThreads, uint32_t tid) {
    (void)tid; (void)nThreads;// to avoid unused warning
    //cerr << "::: tid " << tid << "new" << endl;

    for (;true;) {
        uint64_t ri = gNextIndex++;
        if (ri >= NUM_RELATIONS) return;

        // take the vector with the transactions and sort it by transaction id in order to apply them in order
        auto& relTrans = gTransParseMapPhase[ri];
        if (relTrans.empty()) continue;

        //cerr << "tid " << tid << " got " << ri << " = " << relTrans.size() << endl;

        std::stable_sort(relTrans.begin(), relTrans.end(), TRMapPhaseByTrans);

        auto& relation = gRelations[ri];
        auto& relColumns = gRelColumns[ri].columns;
        uint32_t relCols = gSchema[ri];
        uint64_t lastTransId = relTrans[0].trans_id;
        for (uint32_t c=0; c<relCols; ++c) {
            relColumns[c].transactions.emplace_back(relTrans[0].trans_id, move(vector<CTransStruct>()));
        }
        // for each transaction regarding this relation
        for (auto& trans : relTrans) {
            if (trans.trans_id != lastTransId) {
                // TODO - store the last transactions data
                for (uint32_t c=0; c<relCols; ++c) {
                    sort(relColumns[c].transactions.back().second.begin(), relColumns[c].transactions.back().second.end(), CTransStruct::CompValOnly);
                    // allocate vectors for the current new transaction to put its data
                    relColumns[c].transactions.emplace_back(trans.trans_id, move(vector<CTransStruct>()));
                }
                lastTransId = trans.trans_id;
            }
            if (trans.isDelOp) {
                // this is a delete operation
                vector<tuple_t> operations;
                operations.reserve(trans.rowCount);
                for (const uint64_t* key=trans.values,*keyLimit=key+trans.rowCount;key!=keyLimit;++key) {
                    auto lb = relation.insertedRows.find(*key);
                    if (lb != relation.insertedRows.end()) {
                        // lb->second is a pair<uint64_t, uint64_t*> - trans_id/tuple
                        // - TODO - decrease counter of trans tuples
                        auto tit = lower_bound(relation.transLog.begin(), relation.transLog.end(), lb->second.first, RTLComp);
                        (*tit)->last_del_id = trans.trans_id;
                        --(*tit)->aliveTuples;
                        
                        // update the relation transactions - transfer ownership of the tuple
                        tuple_t tpl = lb->second.second;
                        operations.push_back(tpl);
                        // remove the row from the relations table 
                        relation.insertedRows.erase(lb);
                        for (uint32_t c=0; c<relCols; ++c) {
                            relColumns[c].transactions.back().second.emplace_back(tpl[c], tpl);
                        }
                    }
                }
                if (!operations.empty())
                    relation.transLogDel.emplace_back(trans.trans_id, move(operations));
                delete[] trans.values;
            } else {
                // this is an insert operation
                for (const uint64_t* values=trans.values,*valuesLimit=values+(trans.rowCount*relCols);values!=valuesLimit;values+=relCols) {
                    tuple_t vals = const_cast<uint64_t*>(values);
                    for (uint32_t c=0; c<relCols; ++c) {
                        relColumns[c].transactions.back().second.emplace_back(vals[c], vals);
                    }
                    // finally add the new tuple to the inserted rows of the relation
                    relation.insertedRows[values[0]]=move(make_pair(trans.trans_id, vals));
                }
                // TODO - THIS HAS TO BE IN ORDER - each transaction will have its own transaction history from now on
                relation.transLog.emplace_back(new RelTransLog(trans.trans_id, trans.values, trans.rowCount));
            }
        }
        // TODO - store the last transaction data
        for (uint32_t c=0; c<relCols; ++c) {
            sort(relColumns[c].transactions.back().second.begin(), relColumns[c].transactions.back().second.end(), CTransStruct::CompValOnly);
        }

/*
        // PHASE TWO OF THE ALGORITHM IN THIS STAGE IS TO INCREMENTALLY UPDATE
        // THE INDEXES ONLY FOR THE COLUMNS THAT ARE GOING TO BE REQUESTED IN THE 
        // FOLOWING VALIDATION SESSION - 1st predicates only for now

        // will be used in binary search to find the first column index we need for this relation
        uint32_t rcStart = lp::validation::packRelCol(ri, 0);
        uint32_t rel,col; 
        auto colFrom = lower_bound(gStatColumns->begin(), gStatColumns->end(), rcStart, StatComp);
        if (colFrom == gStatColumns->end()) continue; // no column to index for this relation

        // for each column to be indexed
        for (auto cEnd=gStatColumns->end(); colFrom!=cEnd; ++colFrom) {
            lp::validation::unpackRelCol(colFrom->second, rel, col);
            if (rel != ri) break; // no column to index for this relation
        
            // TODO - IMPORTANT OPTIMIZATION - MERGE THE TWO PASSES INTO ONE

            auto& colTransactions = relColumns[col].transactions;
            auto transFrom = relColumns[col].transTo;
            // first do the inserts - relation.transLog
            for(uint64_t i=0, sz=relation.transLog.size(); i<sz && transFrom<relation.transLog[i]->trans_id; ++i) ++transFrom;
            for(uint64_t tri=0, sz=relation.transLog.size(); tri<sz; ++tri) {
                auto rtl = relation.transLog[tri].get();
                // allocate vectors for the current new transaction to put its data
                colTransactions.emplace_back(rtl->trans_id, move(vector<CTransStruct>()));
                for (uint64_t* values=rtl->tuples.get(),*valuesLimit=values+(rtl->rowCount*relCols);values!=valuesLimit;values+=relCols) {
                    colTransactions.back().second.emplace_back(values[col], values);
                }
            }
            stable_sort(colTransactions.back().second.begin(), colTransactions.back().second.end(), CTransStruct::CompValOnly);

            // second do the deletes - relation.transLogDel
            transFrom = relColumns[col].transTo;
            for(uint64_t i=0, sz=relation.transLogDel.size(); i<sz && transFrom<relation.transLogDel[i].first; ++i) ++transFrom;
            for(uint64_t tri=0, sz=relation.transLogDel.size(); tri<sz; ++tri) {
                auto trp = relation.transLogDel[tri];
                // allocate vectors for the current new transaction to put its data
                colTransactions.emplace_back(trp.first, move(vector<CTransStruct>()));
                for (auto& tpl : trp.second) {
                    colTransactions.back().second.emplace_back(tpl[col], tpl);
                }
            }
            stable_sort(colTransactions.back().second.begin(), colTransactions.back().second.end(), CTransStruct::CompValOnly);
        }
*/

    } // end of while true
}

static inline void checkPendingTransactions(SingleTaskPool& pool) {
#ifdef LPDEBUG
    auto startIndex = LPTimer.getChrono();
#endif

    //cerr << "::: session start ::::" << endl;
    vector<SColType>* cols = &gStats[0].reqCols;
    //for (SColType& cp : *cols) cerr << "==: " << cp.first << " col: " << cp.second << endl; 
    for (uint32_t tid=1; tid<Globals.nThreads; tid++) {
        if (gStats[tid].reqCols.size() > cols->size()) cols = &gStats[tid].reqCols;
    }
    for (uint32_t tid=0; tid<Globals.nThreads; tid++) {
        auto ccols = &gStats[tid].reqCols;
        if (ccols == cols) continue;
        //for (auto& p : *ccols) cerr << " " << p.second;
        copy(ccols->begin(), ccols->end(), back_inserter(*cols));
        ccols->clear();
    }
    std::sort(cols->begin(), cols->end(), StatComp);
    auto it = std::unique(cols->begin(), cols->end(), [] (const SColType& l, const SColType& r) { return l.second == r.second; });
    cols->resize(std::distance(cols->begin(), it));
    //cerr << "unique reqCols: " << cols->size() << endl;
    //for (auto& p : *cols) cerr << " " << p.second;
    //for (SColType& cp : *cols) cerr << "==: " << cp.first << " col: " << cp.second << endl; 
    
    gStatColumns = cols;

    gNextIndex = 0;
    pool.startAll(processPendingIndexTask);
    pool.waitAll();
    for (uint32_t r=0; r<NUM_RELATIONS; ++r) gTransParseMapPhase[r].clear();
    
    // clear the 1st predicate columns 
    cols->clear();

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
    memset(gPendingResults.data(), 0, sizeof(PendingResultType)*gPRsz);
    gNextPending.store(0);

    pool.startAll(processPendingValidationsTask);
    pool.waitAll();

    // update the results - you can get the validation id by adding resIndexOffset to the position
    for (uint64_t i=0, valId=resIndexOffset; i<gPVunique; ++i, ++valId) { 
        gQueryResults.push_back(move(make_pair(valId, gPendingResults[i])));
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
        sort(v.queries.begin(), v.queries.end(), LPQuery::LPQuerySizeLessThan);

        bool conflict = false, otherFinishedThis = false;
        // for each query in this validation         
        for (auto& q : v.queries) {
            if (atoRes) { otherFinishedThis = true; /*cerr << "h" << endl;*/ break; }

            auto& relColumns = gRelColumns[q.relationId].columns;

            
            // protect from the case where there is no single predicate
            if (q.predicates.empty()) { 
                //cerr << "empty: " << v.validationId << endl; 
                auto& transactionsCheck = relColumns[0].transactions;
                auto transFromCheck = std::lower_bound(transactionsCheck.begin(), transactionsCheck.end(), v.from, CTRSLessThan);
                if (transFromCheck == transactionsCheck.end()) { 
                    // no transactions exist for this query
                    continue;
                } else { 
                    // transactions exist for this query
                    conflict=true;  break; 
                }; 
            }

            // IMPORTANT!!! - sort them in order to have the equality checks first
            std::sort(q.predicates.begin(), q.predicates.end(), LPQuery::QCSortOp);

            auto& pFirst = q.predicates[0];
            auto& transactions = relColumns[pFirst.column].transactions;
            auto transFrom = std::lower_bound(transactions.begin(), transactions.end(), v.from, CTRSLessThan);
            auto transTo = std::upper_bound(transactions.begin(), transactions.end(), v.to, CTRSLessThan);

            //cerr << "after: " << v.from << "-" << v.to << "=" << (transTo-transFrom) << " for col: " << pFirst.column << "-" << pFirst.value << endl;

            uint32_t pFrom = 0;
            for(auto iter=transFrom; iter!=transTo; ++iter) {  
                auto& transValues = iter->second;
                decltype(transValues.begin()) tupFrom, tupTo, 
                    tBegin = transValues.begin(), tEnd=transValues.end();
                // find the valid tuples using range binary searches based on the first predicate
                if (pFirst.op == lp::LPOps::Equal) {
                    tupFrom = std::lower_bound(tBegin, tEnd, pFirst.value, CTRSValueLessThan);
                    if (tupFrom == tEnd) continue;
                    tupTo = std::upper_bound(tBegin, tEnd, pFirst.value, CTRSValueLessThan);                   
                    pFrom = 1;
                } else if (pFirst.op == lp::LPOps::Less) {
                    tupFrom = tBegin;                    
                    tupTo = std::lower_bound(tBegin, tEnd, pFirst.value, CTRSValueLessThan);                   
                    pFrom = 1;
                } else if (pFirst.op == lp::LPOps::LessOrEqual) {
                    tupFrom = tBegin;                    
                    tupTo = std::upper_bound(tBegin, tEnd, pFirst.value, CTRSValueLessThan);                   
                    pFrom = 1;
                } else if (pFirst.op == lp::LPOps::Greater) {
                    tupFrom = std::upper_bound(tBegin, tEnd, pFirst.value, CTRSValueLessThan);  
                    tupTo = tEnd;                   
                    pFrom = 1;
                } else if (pFirst.op == lp::LPOps::GreaterOrEqual) {
                    tupFrom = std::lower_bound(tBegin, tEnd, pFirst.value, CTRSValueLessThan);
                    tupTo = tEnd;                   
                    pFrom = 1;
                } else {
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
                        //if (v.validationId == 4) cerr << "pred:" << c << endl;
                        // make the actual check
                        uint64_t tupleValue = tuple[c.column]; 
                        //cerr << "tpl value: " << tupleValue << endl;
                        uint64_t queryValue = c.value;
                        bool result=false;
                        switch (c.op) {
                            case lp::LPOps::Equal: 
                                result=(tupleValue==queryValue); 
                                break;
                            case lp::LPOps::Less: 
                                result=(tupleValue<queryValue); 
                                break;
                            case lp::LPOps::LessOrEqual: 
                                result=(tupleValue<=queryValue); 
                                break;
                            case lp::LPOps::Greater: 
                                result=(tupleValue>queryValue); 
                                break;
                            case lp::LPOps::GreaterOrEqual: 
                                result=(tupleValue>=queryValue); 
                                break;
                            case lp::LPOps::NotEqual: 
                            case lp::LPOps::NotEqualLast: 
                                result=(tupleValue!=queryValue); 
                                break;
                        } 
                        // there is one predicate not true so this whole query on this relation is false
                        if (!result) { match=false; break; }
                    } // end of single query predicates
                    //cerr << "match: " << match << " conflict: " << conflict << endl;
                    if (match) { conflict=true; break; }    
                } // end of all tuples for this transaction
                if (conflict) { break; }
            } // end of all the transactions for this relation for this specific query
            if (conflict) break;
        }// end for all queries
        // update the pending results to help other threads skip this validation 
        // if it has other parts
        if (conflict && !otherFinishedThis)  { /*cerr<< "c: " << v.validationId << endl;*/  atoRes = true;}
    } // while true take more validations 
}

