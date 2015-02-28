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

#include "ReferenceTypes.hpp"
#include "LPQueryTypes.hpp"

#include "atomicwrapper.hpp"
#include "LPTimer.hpp"
#include "BoundedQueue.hpp"
#include "BoundedAlloc.hpp"
#include "SingleTaskPool.hpp"
#include "MultiTaskPool.hpp"


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
    os << "{" /*<< o.trans_id*/ << "-" << o.value << "}";
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
    //vector<CTransStruct> transactions;
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

// the structure that WAS used inside the RelationStruct to record the transactoins of that relation
struct TransStruct {
    uint64_t trans_id;
    std::unique_ptr<tuple_t> tuple;

    TransStruct(uint64_t trid, std::unique_ptr<tuple_t> t):
        trans_id(trid), tuple(move(t)) {}
};
struct TRSLessThan_t {
    bool operator() (const TransStruct& left, const TransStruct& right) {
        return left.trans_id < right.trans_id;
    }
    bool operator() (const TransStruct& o, uint64_t target) {
        return o.trans_id < target;
    }
    bool operator() (uint64_t target, const TransStruct& o) {
        return target < o.trans_id;
    }
} TRSLessThan;

// The general structure for each relation
struct RelationStruct {
    unordered_map<uint32_t, std::unique_ptr<uint64_t[]>> insertedRows;
};

// RELATION STRUCTURES

// general relations
static std::unique_ptr<RelationStruct[]> gRelations;
// transactions in each relation column - all tuples of same transaction in one vector
struct RelationColumns {
    vector<ColumnStruct> columns;
};
static std::unique_ptr<RelationColumns[]> gRelColumns;



struct TRMapPhase {
    uint64_t trans_id;
    bool isDelOp;
    uint32_t rowCount;
    uint64_t *values; // delete op => row keys to delete | insert => tuples

    TRMapPhase(uint64_t tid, bool isdel, uint32_t rows, uint64_t *vals)
        : trans_id(tid), isDelOp(isdel), rowCount(rows), values(vals) {}
};
static unique_ptr<std::mutex[]> gRelTransMutex;
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
    vector<TransOperation> operations;
    TransactionStruct(uint64_t tid, vector<TransOperation> ops) : trans_id(tid), operations(move(ops)) {}  
};
static vector<TransactionStruct> gTransactionHistory;
static std::mutex gTransactionHistoryMutex;

static uint32_t NUM_RELATIONS;
static std::unique_ptr<uint32_t[]> gSchema;

///////// AUXILIARY STRUCTURES FOR THE WHOLE PROGRAM

static vector<pair<uint64_t, unique_ptr<map<uint32_t, vector<tuple_t>>>>> gPendingIndex;
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
struct ParseValidationStruct {
    ReceivedMessage *msg;
    uint64_t refId;
    BoundedQueue<ReceivedMessage> *msgQ;
    uint64_t memRefId;
    BoundedAlloc<ParseValidationStruct> *memQ;
};


LPTimer_t LPTimer;

struct GlobalState {
    enum State : uint32_t { SCHEMA, TRANSACTION, VALIDATION, FORGET, FLUSH };
    State state;
} Globals;

//---------------------------------------------------------------------------

// JUST SOME FUNCTION DECLARATIONS THAT ARE DEFINED BELOW
static void checkPendingValidations(SingleTaskPool&);
static void processPendingValidationsTask(uint32_t nThreads, uint32_t tid);
static void processSingleTransaction(const Transaction&);

static inline void checkPendingTransactions(SingleTaskPool& pool);
static void processPendingIndexTask(uint32_t nThreads, uint32_t tid);

static void processTransactionMessage(const Transaction& t, vector<char>& data);
///--------------------------------------------------------------------------
///--------------------------------------------------------------------------

static void processDefineSchema(const DefineSchema& d) {
    gSchema.reset(new uint32_t[d.relationCount]);
    memcpy(gSchema.get(), d.columnCounts, sizeof(uint32_t)*d.relationCount);
    NUM_RELATIONS = d.relationCount;

    gRelTransMutex.reset(new std::mutex[d.relationCount]);
    gTransParseMapPhase.reset(new vector<TRMapPhase>[d.relationCount]);

    gRelations.reset(new RelationStruct[d.relationCount]);
    gRelColumns.reset(new RelationColumns[d.relationCount]);
    for(uint32_t ci=0; ci<d.relationCount; ++ci) {
        gRelColumns[ci].columns.resize(gSchema[ci]);
    }
}
//---------------------------------------------------------------------------

static inline void processTransaction(const Transaction& t, const vector<char>& tdata) {
    (void)tdata;
    processSingleTransaction(t);
    // hack to get the pointer to the beginning of the Transactiob structure!!!
    //cerr << "tr:" << t.transactionId << endl;
    //gPendingTransactions.push_back(reinterpret_cast<const Transaction*>((char*)t.operations-sizeof(uint64_t)-2*sizeof(uint32_t)));
    //gPendingTransactions.push_back(t);
}

//---------------------------------------------------------------------------


static void processValidationQueries(const ValidationQueries& v, const vector<char>& vdata) {
#ifdef LPDEBUG
    auto start = LPTimer.getChrono();    
#endif

    (void)vdata;
    // TODO - OPTIMIZATION CAN BE DONE IF I JUST COPY THE WHOLE DATA instead of parsing it
    // try to put all the queries into a vector
    vector<LPQuery> queries;
    const char* qreader=v.queries;
    const Query *q;
    for (unsigned int i=0;i<v.queryCount;++i) {
        q=reinterpret_cast<const Query*>(qreader);
        LPQuery nQ(*q);
        //cerr << v.validationId << "====" << v.from << ":" << v.to << nQ << endl;
        if (!lp::validation::isQueryUnsolvable(nQ)) queries.push_back(move(nQ));
        //queries.push_back(move(LPQuery(*q)));
        qreader+=sizeof(Query)+(sizeof(Query::Column)*q->columnCount);
    }
    //if (v.validationId == 19631)
    //  cerr << v.validationId << "====" << v.from << ":" << v.to << "=" << v.queryCount << "=" << queries << endl;

    uint64_t batchPos = v.from; 

    //uint64_t trRange = v.to - v.from + 1;
    //static uint64_t bSize = 500;
    {
        std::lock_guard<std::mutex> lk(gPendingValidationsMutex);

        //  while (trRange > bSize) {
        //cerr << batchPos << "-" << batchPos+bSize-1 << endl;
        //gPendingValidations.push_back(move(LPValidation(v.validationId, batchPos, batchPos+bSize-1, queries)));    
        //trRange -= bSize; batchPos += bSize;
        //}
        gPendingValidations.push_back(move(LPValidation(v.validationId, batchPos, v.to, move(queries))));    
        //cerr << batchPos << "-" << v.to << endl;

        // update the global pending validations to reflect this new one
        ++gPVunique;
    }
    //cerr << "Success Validate: " << v.validationId << " ::" << v.from << ":" << v.to << " result: " << conflict << endl;
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
        for (auto& cCol : cRelCol.columns) {
            cCol.transactions.erase(cCol.transactions.begin(),
                    upper_bound(cCol.transactions.begin(), cCol.transactions.end(),
                        f.transactionId,
                        [](const uint64_t target, const ColumnTransaction_t& ct){ return target < ct.first; }
                        ));
        }
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
    ParseValidationStruct *pvs = static_cast<ParseValidationStruct*>(args);
    processValidationQueries(*reinterpret_cast<const ValidationQueries*>(pvs->msg->data.data()), pvs->msg->data); 
    pvs->msgQ->registerDeq(pvs->refId);
    pvs->memQ->free(pvs->memRefId);
}


void inline parseTransactionPH1(uint32_t nThreads, uint32_t tid, void *args) {
    (void)tid; (void)nThreads;
    ParseValidationStruct *pvs = static_cast<ParseValidationStruct*>(args);
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

    uint64_t MessageQSize = 500;
    BoundedQueue<ReceivedMessage> msgQ(MessageQSize);
    BoundedAlloc<ParseValidationStruct> memQ(MessageQSize);

    std::thread readerTask(ReaderTask, std::ref(msgQ));

    SingleTaskPool workerThreads(numOfThreads, processPendingValidationsTask);
    workerThreads.initThreads();

    // leave two available workes - master - msgQ
    MultiTaskPool multiPool(1);
    multiPool.initThreads();
    multiPool.startAll();

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
#ifndef LPDEBUG
                        ++gTotalValidations; // this is just to count the total validations....not really needed!
#endif
                        //ParseValidationStruct *pvs = new ParseValidationStruct();
                        BoundedAlloc<ParseValidationStruct>::BAResult& mem = memQ.malloc();
                        ParseValidationStruct *pvs = mem.value;
                        pvs->msgQ = &msgQ;
                        pvs->refId = res.refId;
                        pvs->memRefId = mem.refId;
                        pvs->memQ = &memQ;
                        pvs->msg = &msg;
                        multiPool.addTask(parseValidation, static_cast<void*>(pvs)); 
                        break;
                    }
                case MessageHead::Transaction: 
                    {Globals.state = GlobalState::TRANSACTION;
                    //processTransaction(*reinterpret_cast<const Transaction*>(msg.data.data()), msg.data); 
                    //msgQ.registerDeq(res.refId);
                    BoundedAlloc<ParseValidationStruct>::BAResult& mem = memQ.malloc();
                    ParseValidationStruct *pvs = mem.value;
                    pvs->msgQ = &msgQ;
                    pvs->refId = res.refId;
                    pvs->memRefId = mem.refId;
                    pvs->memQ = &memQ;
                    pvs->msg = &msg;
                    multiPool.addTask(parseTransactionPH1, static_cast<void*>(pvs)); 
                    break;
                    }
                case MessageHead::Flush:  
                    // check if we have pending transactions to be processed
                    multiPool.helpExecution();
                    multiPool.waitAll();
                    checkPendingValidations(workerThreads);
                    Globals.state = GlobalState::FLUSH;
                    processFlush(*reinterpret_cast<const Flush*>(msg.data.data())); 
                    msgQ.registerDeq(res.refId);
                    break;

                case MessageHead::Forget: 
                    // check if we have pending transactions to be processed
                    multiPool.helpExecution();
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


static void processSingleTransaction(const Transaction& t) {
#ifdef LPDEBUG
    auto start = LPTimer.getChrono();
    ++gTotalTransactions; 
#endif
    const char* reader=t.operations;

    vector<TransOperation> operations;
    unique_ptr<map<uint32_t, vector<tuple_t>>> locals_ptr(new map<uint32_t, vector<tuple_t>>);
    auto& locals = *locals_ptr;

    // Delete all indicated tuples
    for (uint32_t index=0;index!=t.deleteCount;++index) {
        auto& o=*reinterpret_cast<const TransactionOperationDelete*>(reader);
        auto& rows = gRelations[o.relationId].insertedRows;

        // TODO - lock here to make it to make all the deletions parallel naive locking first - 
        // TODO try to lock with try_lock and try again at the end if some relations failed
        {// start of lock_guard
            //std::lock_guard<std::mutex> lk(gRelTransMutex[o.relationId]);
            auto& vec = locals[o.relationId];
            for (const uint64_t* key=o.keys,*keyLimit=key+o.rowCount;key!=keyLimit;++key) {
                auto lb = rows.find(*key);
                if (lb != rows.end()) {
                    // update the relation transactions - transfer ownership of the tuple
                    operations.push_back(move(TransOperation(o.relationId, move(lb->second))));
                    vec.push_back(operations.back().tuple.get());

                    // remove the row from the relations table
                    rows.erase(lb);
#ifdef LPDEBUG
                    ++gTotalTuples;
#endif
                }
            }
        }// end of lock_guard

        // advance to the next Relation deletions
        reader+=sizeof(TransactionOperationDelete)+(sizeof(uint64_t)*o.rowCount);
    }

    // Insert new tuples
    for (uint32_t index=0;index!=t.insertCount;++index) {
        auto& o=*reinterpret_cast<const TransactionOperationInsert*>(reader);
        auto& relation = gRelations[o.relationId];
        const uint32_t relCols = gSchema[o.relationId];
        // TODO - lock here to make it to make all the deletions parallel naive locking first - 
        // TODO try to lock with try_lock and try again at the end if some relations failed
        {// start of lock_guard
            //std::lock_guard<std::mutex> lk(gRelTransMutex[o.relationId]);
            auto& vec = locals[o.relationId];
            uint64_t* tptr_r;
            const uint64_t *vptr;
            for (const uint64_t* values=o.values,*valuesLimit=values+(o.rowCount*relCols);values!=valuesLimit;values+=relCols) {
                std::unique_ptr<uint64_t[]> tptr2(new uint64_t[relCols]);
                tptr_r = tptr2.get(); vptr = values;
                for (uint32_t c=0; c<relCols; ++c) *tptr_r++ = *vptr++;

                // store the pointer to the tuple to the transaction actions locally
                vec.push_back(tptr2.get());

                // finally add the new tuple to the inserted rows of the relation
                relation.insertedRows[values[0]]=move(tptr2);

#ifdef LPDEBUG
                ++gTotalTuples;
#endif
            }
        }// end of lock_guard
        // advance to next Relation insertions
        reader+=sizeof(TransactionOperationInsert)+(sizeof(uint64_t)*o.rowCount*relCols);
    }

    gPendingIndex.push_back(move(make_pair(t.transactionId, move(locals_ptr))));

    // update the transaction history
    // TODO - HERE WE WILL HAVE TO LOCK THE VECTOR AND ADD THE TRANSACTION IN THE RIGHT PLACE
    gTransactionHistory.push_back(move(TransactionStruct(t.transactionId, move(operations))));


#ifdef LPDEBUG
    LPTimer.transactions += LPTimer.getChrono(start);
#endif
}


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
    ++gTotalTransactions; 
#endif
    (void)data;

    const char* reader=t.operations;
    // Delete all indicated tuples
    for (uint32_t index=0;index!=t.deleteCount;++index) {
        auto& o=*reinterpret_cast<const TransactionOperationDelete*>(reader);
        // TODO - lock here to make it to make all the deletions parallel naive locking first - 
        // TODO try to lock with try_lock and try again at the end if some relations failed
        {// start of lock_guard
            std::lock_guard<std::mutex> lk(gRelTransMutex[o.relationId]);
            uint64_t *ptr = new uint64_t[o.rowCount];
            const uint64_t *keys = o.keys;
            for (uint32_t c=0; c<o.rowCount; ++c) *ptr++ = *keys++;
            gTransParseMapPhase[o.relationId].emplace_back(t.transactionId, true, o.rowCount, ptr);
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
            std::lock_guard<std::mutex> lk(gRelTransMutex[o.relationId]);
            uint64_t* tptr = new uint64_t[relCols*o.rowCount];
            const uint64_t *vptr=o.values;
            uint32_t sz = relCols*o.rowCount;
            for (uint32_t c=0; c<sz; ++c) *tptr++ = *vptr++;
            gTransParseMapPhase[o.relationId].emplace_back(t.transactionId, false, o.rowCount, tptr-sz);
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
    return l.trans_id < r.trans_id;
}

static void processPendingIndexTask(uint32_t nThreads, uint32_t tid) {
    (void)tid; (void)nThreads;// to avoid unused warning
    for (;true;) {
        uint64_t ri = gNextIndex++;
        if (ri >= NUM_RELATIONS) return;

        // take the vector with the transactions and sort it by transaction id in order to apply them in order
        auto& relTrans = gTransParseMapPhase[ri];
        if (relTrans.empty()) continue;

        cerr << "tid " << tid << " got " << ri << " = " << relTrans.size() << endl;

        std::sort(relTrans.begin(), relTrans.end(), TRMapPhaseByTrans);

        auto& relation = gRelations[ri];
        auto& relColumns = gRelColumns[ri];
        uint32_t relCols = gSchema[ri];
        uint64_t lastTransId = relTrans[0].trans_id;
        std::unique_ptr<vector<CTransStruct>[]> colPtrs(new vector<CTransStruct>[relCols]);
        for (uint32_t c=0; c<relCols; ++c) {
            relColumns.columns[c].transactions.emplace_back(move(make_pair(relTrans[0].trans_id, move(vector<CTransStruct>()))));
        }
        // for each transaction regarding this relation
        for (auto& trans : relTrans) {
            if (trans.trans_id != lastTransId) {
                // TODO - store the last transactions data
                for (uint32_t c=0; c<relCols; ++c) {
                    sort(relColumns.columns[c].transactions.back().second.begin(), relColumns.columns[c].transactions.back().second.end(), CTransStruct::CompValOnly);
                    // allocate vectors for the current new transaction to put its data
                    relColumns.columns[c].transactions.emplace_back(move(make_pair(trans.trans_id, move(vector<CTransStruct>()))));
                }
            }
            if (trans.isDelOp) {
                // this is a delete operation
                vector<TransOperation> operations;
                for (const uint64_t* key=trans.values,*keyLimit=key+trans.rowCount;key!=keyLimit;++key) {
                    auto lb = relation.insertedRows.find(*key);
                    if (lb != relation.insertedRows.end()) {
                        // update the relation transactions - transfer ownership of the tuple
                        operations.push_back(move(TransOperation(ri, move(lb->second))));
                        // remove the row from the relations table
                        relation.insertedRows.erase(lb);
                        for (uint32_t c=0; c<relCols; ++c) {
                            relColumns.columns[c].transactions.back().second.emplace_back(operations.back().tuple[c], operations.back().tuple.get());
                        }
                    }
                }
                {
                    std::lock_guard<std::mutex> lk(gTransactionHistoryMutex);
                    gTransactionHistory.push_back(move(TransactionStruct(trans.trans_id, move(operations))));
                }
            } else {
                // this is an insert operation
                uint64_t *tptr_r; const uint64_t *vptr;
                for (const uint64_t* values=trans.values,*valuesLimit=values+(trans.rowCount*relCols);values!=valuesLimit;values+=relCols) {
                    std::unique_ptr<uint64_t[]> tptr2(new uint64_t[relCols]);
                    tptr_r = tptr2.get(); vptr = values;
                    for (uint32_t c=0; c<relCols; ++c) {
                        *tptr_r++ = *vptr++;
                        relColumns.columns[c].transactions.back().second.emplace_back(values[c], tptr2.get());
                        cerr << values[c] << "-" << tptr2[c] << ":" << relColumns.columns[c].transactions.back().second.back().value << endl;
                    }
                    // finally add the new tuple to the inserted rows of the relation
                    relation.insertedRows[values[0]]=move(tptr2);
                } 
            }
            delete[] trans.values;
        }
        // TODO - store the last transaction data
        for (uint32_t c=0; c<relCols; ++c) {
            sort(relColumns.columns[c].transactions.back().second.begin(), relColumns.columns[c].transactions.back().second.end(), CTransStruct::CompValOnly);
        }


/*
        ///// MOST DIFFICULT PART - NEED TO BE OPTIMIZED - Insert the tuples in the relations
        // column-wise !!!
        auto& trans_locals = gPendingIndex[ii];
        for (auto& rtr : *trans_locals.second) {
            uint32_t relCols = gSchema[rtr.first];
            std::unique_ptr<vector<CTransStruct>[]> colPtrs(new vector<CTransStruct>[relCols]);
            for (uint32_t ti=0,tsz=rtr.second.size(); ti<tsz; ++ti) {
                for (uint32_t c=0; c<relCols; ++c) {
                    colPtrs[c].push_back(move(CTransStruct(rtr.second[ti][c], rtr.second[ti])));
                }
            }
            for (uint32_t c=0; c<relCols; ++c) {
                sort(colPtrs[c].begin(), colPtrs[c].end(), CTransStruct::CompValOnly);
                //gRelColumns[rtr.first].columns[c].transactions.push_back(move(std::make_pair(t.transactionId, move(colPtrs[c]))));
            }

            // INSERT the vectors into the relation columns in the right position to retain transaction ordering
            {// start of lock_guard
                std::lock_guard<std::mutex> lk(gRelTransMutex[rtr.first]);
                auto& trans = gRelColumns[rtr.first].columns[0].transactions;
                auto lb = lower_bound(trans.begin(), trans.end(), trans_locals.first, 
                        [](const pair<uint64_t, vector<CTransStruct>>& o, uint64_t target) { return o.first < target; }
                        );
                uint64_t diff = lb - trans.begin();
                // TODO - find the position that we need to get inserted - BINARY SEARCH IF PARALLEL
                for (uint32_t c=0; c<relCols; ++c) {
                    //gRelColumns[rtr.first].columns[c].transactions.push_back(move(std::make_pair(trans_locals.first, move(colPtrs[c]))));
                    auto& tr = gRelColumns[rtr.first].columns[c].transactions;
                    tr.insert(tr.begin()+diff, move(std::make_pair(trans_locals.first, move(colPtrs[c]))));
                }
            }

        }
*/
    } // end of while true
}

static inline void checkPendingTransactions(SingleTaskPool& pool) {
#ifdef LPDEBUG
    auto startIndex = LPTimer.getChrono();
#endif
    gNextIndex = 0;
    pool.startAll(processPendingIndexTask);
    pool.waitAll();
    for (uint32_t r=0; r<NUM_RELATIONS; ++r) gTransParseMapPhase[r].clear();
#ifdef LPDEBUG
    LPTimer.transactionsIndex += LPTimer.getChrono(startIndex);
#endif
}


static uint64_t resIndexOffset = 0;
static std::atomic<uint64_t> gNextPending;

static void checkPendingValidations(SingleTaskPool &pool) {
    if (gPendingValidations.empty()) return;
#ifdef LPDEBUG
    auto start = LPTimer.getChrono();
#endif
    checkPendingTransactions(pool);

    //cerr << gPendingValidations.size() << " " << endl;
    // find the min & max validation id
    // assuming that gPendingValidations is sorted on the validation Id
    resIndexOffset = UINT64_MAX;
    for (auto& pv : gPendingValidations) if (pv.validationId < resIndexOffset) resIndexOffset = pv.validationId;
    auto gPRsz = gPendingResults.size();
    if (gPVunique > gPRsz)
        gPendingResults.resize(gPVunique);
    memset(gPendingResults.data(), 0, sizeof(PendingResultType)*gPRsz);
    gNextPending.store(0);

    //cerr << "before start!" << endl;
    pool.startAll(processPendingValidationsTask);
    //processPendingValidations();
    //cerr << "before wait!" << endl;
    pool.waitAll();
    //cerr << "after wait!" << endl;

    // update the results - you can get the validation id by adding resIndexOffset to the position
    for (uint64_t i=0, valId=resIndexOffset; i<gPVunique; ++i, ++valId) { 
        //cerr << valId << ":" << gPendingResults[valId-resIndexOffset].load() << " ";
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
        // check if the query is by default false - NON-CONFLICT
        if (lp::validation::unsolvable(v)) {
        //cerr << "unsolvable" << endl;
        continue;
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

