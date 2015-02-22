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

#include "LPTimer.hpp"
//#include "BoundedSRSWQueue.hpp"
#include "BoundedQueue.hpp"
#include "SingleTaskPool.hpp"

//---------------------------------------------------------------------------
using namespace std;

#define LPDEBUG  

//---------------------------------------------------------------------------
// Wire protocol messages
//---------------------------------------------------------------------------
struct MessageHead {
    /// Message types
    enum Type : uint32_t { Done, DefineSchema, Transaction, ValidationQueries, Flush, Forget };
    /// Total message length, excluding this head
    uint32_t messageLen;
    /// The message type
    Type type;
};
//---------------------------------------------------------------------------
struct DefineSchema {
    /// Number of relations
    uint32_t relationCount;
    /// Column counts per relation, one count per relation. The first column is always the primary key
    uint32_t columnCounts[];
};
//---------------------------------------------------------------------------
struct Transaction {
    /// The transaction id. Monotonic increasing
    uint64_t transactionId;
    /// The operation counts
    uint32_t deleteCount,insertCount;
    /// A sequence of transaction operations. Deletes first, total deleteCount+insertCount operations
    char operations[];
};
//---------------------------------------------------------------------------
struct TransactionOperationDelete {
    /// The affected relation
    uint32_t relationId;
    /// The row count
    uint32_t rowCount;
    /// The deleted values, rowCount primary keyss
    uint64_t keys[];
};
//---------------------------------------------------------------------------
struct TransactionOperationInsert {
    /// The affected relation
    uint32_t relationId;
    /// The row count
    uint32_t rowCount;
    /// The inserted values, rowCount*relation[relationId].columnCount values
    uint64_t values[];
};
//---------------------------------------------------------------------------
struct ValidationQueries {
    /// The validation id. Monotonic increasing
    uint64_t validationId;
    /// The transaction range
    uint64_t from,to;
    /// The query count
    uint32_t queryCount;
    /// The queries
    char queries[];
};
//---------------------------------------------------------------------------
struct Query {
    /// A column description
    struct Column {
        /// Support operations
        enum Op : uint32_t { Equal, NotEqual, Less, LessOrEqual, Greater, GreaterOrEqual };
        /// The column id
        uint32_t column;
        /// The operations
        Op op;
        /// The constant
        uint64_t value;
    };

    /// The relation
    uint32_t relationId;
    /// The number of bound columns
    uint32_t columnCount;
    /// The bindings
    Column columns[];
};
//---------------------------------------------------------------------------
struct Flush {
    /// All validations to this id (including) must be answered
    uint64_t validationId;
};
//---------------------------------------------------------------------------
struct Forget {
    /// Transactions older than that (including) will not be tested for
    uint64_t transactionId;
};
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

struct LPQuery {
    /// The relation
    uint32_t relationId;
    /// The number of bound columns
    uint32_t columnCount;
    /// The bindings
    vector<Query::Column> predicates;
    LPQuery() : relationId(-1), columnCount(0), predicates(vector<Query::Column>()) {}
    LPQuery(const Query& q) : relationId(q.relationId), columnCount(q.columnCount), predicates(vector<Query::Column>(q.columnCount)) {
        memcpy(predicates.data(), q.columns, sizeof(Query::Column)*columnCount);
        std::sort(predicates.begin(), predicates.end());
        predicates.resize(std::distance(predicates.begin(), std::unique(predicates.begin(), predicates.end())));
        //if (columns.size() != columnCount) cerr << "diff: " << columnCount-columns.size() << endl;
        // sort them in order to have the equality checks first
        std::sort(predicates.begin(), predicates.end(), LPQuery::QCSortOp);
        columnCount = predicates.size();
    }
    static bool QCSortOp (const Query::Column& left, const Query::Column& right) {
        if (left.op < right.op) return true;
        else if (right.op < left.op) return false;
        else if (left.column < right.column) return true;
        else if (right.column < left.column) return false;
        else return left.value < right.value;    
    }
    static bool LPQuerySizeLessThan(const LPQuery& left, const LPQuery& right) {
        if (left.columnCount < right.columnCount) return true;
        else if (right.columnCount < left.columnCount) return false;
        else return left.relationId < right.relationId;
    }
};
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
ostream& operator<< (ostream& os, const LPQuery& o) {
    os << "{" << o.relationId << "-" << o.columnCount << ":: " << o.predicates << "}";
    return os;
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
ostream& operator<< (ostream& os, const Query::Column& o) {
    os << "[" << o.column << ":" << o.op << ":" << o.value << "]";
    return os;
}

// Custom data structures to hold data
struct ColumnStruct {
    //vector<CTransStruct> transactions;
    uint64_t minVal;
    uint64_t maxVal;
};
struct TransStruct {
    uint64_t trans_id;
    vector<uint64_t> tuple;
    TransStruct(uint64_t trid, vector<uint64_t> t):
        trans_id(trid), tuple(t) {}
    TransStruct(uint64_t trid):
        trans_id(trid), tuple(vector<uint64_t>()) {}
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
struct RelationStruct {
    vector<ColumnStruct> columns;
    vector<TransStruct> transactions;
    unordered_map<uint32_t, vector<uint64_t>> insertedRows;
};
static vector<RelationStruct> gRelations;

struct TransOperation {
    uint32_t rel_id;
    uint64_t row_id;
    TransOperation(uint32_t relid, uint64_t rowid):
        rel_id(relid), row_id(rowid) {}
};
struct TransactionStruct {
    //    uint64_t trans_id; // we will use direct indexing
    vector<TransOperation> operations;
    TransactionStruct(vector<TransOperation> ops) : operations(ops) {}
};

static vector<uint32_t> gSchema;


template <typename T>
class atomic_wrapper {
    private:
        std::atomic<T> _a;
        //char _padding[64-sizeof(std::atomic<T>)];
    public:
        atomic_wrapper()
            :_a() {}

        atomic_wrapper(const std::atomic<T> &a)
            :_a(a.load()){}

        atomic_wrapper(const atomic_wrapper &other)
            :_a(other._a.load()){}

        atomic_wrapper &operator=(const atomic_wrapper &other) {
            _a.store(other._a.load());
            return *this;
        }

        void store(T desired) {_a.store(desired);}
        void store(T desired) volatile {_a.store(desired);}
        T load() const { return _a.load(); }
        T load() const volatile { return _a.load(); }

        std::atomic<T>& operator=(T desired) { _a.store(desired); return _a; }
        std::atomic<T>& operator=(T desired) volatile { _a.store(desired); return _a; }

        operator bool() const { return _a.load(); }
};
struct LPValidation {
    uint64_t validationId;
    uint64_t from,to;
    vector<LPQuery> queries;
    LPValidation(const ValidationQueries& v, vector<LPQuery> q)
        : validationId(v.validationId), from(v.from), to(v.to), queries(move(q)) {}
    LPValidation(uint64_t vid, uint64_t fr, uint64_t t, vector<LPQuery> q)
        : validationId(vid), from(fr), to(t), queries(q) {}
};
static vector<LPValidation> gPendingValidations;
//typedef atomic_wrapper<bool> PendingResultType;
typedef char PendingResultType;
static vector<PendingResultType> gPendingResults;
static vector<pair<uint64_t,bool>> gQueryResults;
static uint64_t gPVunique;

template<class T> using SRSWQueue = BoundedQueue<T>;
struct ReceivedMessage {
    MessageHead head;
    vector<char> data;
};
static vector<SRSWQueue<ReceivedMessage>::BQResult> gPendingTransactions;

static std::mutex *gRelTransMutex;


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

static uint64_t processPendingMessages(SingleTaskPool&, SRSWQueue<ReceivedMessage>&);
static uint64_t processPendingTransactions(SingleTaskPool&);
static void processPendingTransactionsTask(uint32_t nThreads, uint32_t tid);



///--------------------------------------------------------------------------
///--------------------------------------------------------------------------

static void processDefineSchema(const DefineSchema& d) {
    gSchema.clear();
    gSchema.insert(gSchema.begin(),d.columnCounts,d.columnCounts+d.relationCount);

    gRelations.clear();
    gRelations.resize(d.relationCount);
    for(uint32_t ci=0; ci<d.relationCount; ++ci) {
        gRelations[ci].columns.resize(gSchema[ci]);
    }
}
//---------------------------------------------------------------------------
/*
static void processTransaction(vector<char>& data) {
    gPendingTransactions.push_back(data.data());
}
*/
/*
static void processTransaction(const Transaction& t) {
//static void processTransaction(char* t) {
    processSingleTransaction(t);
    // hack to get the pointer to the beginning of the Transactiob structure!!!
    //cerr << "tr:" << t.transactionId << endl;
    //gPendingTransactions.push_back(reinterpret_cast<const Transaction*>((char*)t.operations-sizeof(uint64_t)-2*sizeof(uint32_t)));
    //gPendingTransactions.push_back(t);
}
*/

//---------------------------------------------------------------------------
static void processValidationQueries(const ValidationQueries& v) {
#ifdef LPDEBUG
    auto start = LPTimer.getChrono();    
#endif

    // TODO - OPTIMIZATION CAN BE DONE IF I JUST COPY THE WHOLE DATA instead of parsing it
    // try to put all the queries into a vector
    vector<LPQuery> queries;
    const char* qreader=v.queries;
    const Query *q;
    for (unsigned int i=0;i<v.queryCount;++i) {
        q=reinterpret_cast<const Query*>(qreader);
        queries.push_back(move(LPQuery(*q)));
        qreader+=sizeof(Query)+(sizeof(Query::Column)*q->columnCount);
    }
    // sort the queries based on everything to remove duplicates
    //sort(queries.begin(), queries.end());
    //cerr << "size before: " << queries.size() << endl;
    //queries.resize(std::distance(queries.begin(), unique(queries.begin(), queries.end())));
    //cerr << "size after: " << queries.size() << endl;
    // sort the queries based on the number of the columns needed to check
    // small queries first in order to try finding a solution faster
    sort(queries.begin(), queries.end(), LPQuery::LPQuerySizeLessThan);
    //cerr << "====" << v.from << ":" << v.to << endl;
    uint64_t batchPos = v.from; 
    
    uint64_t trRange = v.to - v.from + 1;
    uint64_t bSize = 1000;
    while (trRange > bSize) {
        //cerr << batchPos << "-" << batchPos+bSize-1 << endl;
        gPendingValidations.push_back(move(LPValidation(v.validationId, batchPos, batchPos+bSize-1, queries)));    
        trRange -= bSize; batchPos += bSize;
    }
    
    gPendingValidations.push_back(move(LPValidation(v.validationId, batchPos, v.to, move(queries))));    
    //cerr << batchPos << "-" << v.to << endl;

    // update the global pending validations to reflect this new one
    ++gPVunique;

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
    char zero = 48;
    char one = 49;
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

    TransStruct fstruct(f.transactionId);
    // Delete the transactions from inside the columns in the relations
    for(auto crel=gRelations.begin(), iend=gRelations.end(); crel!=iend; ++crel) {
        // delete this transaction from the lastRel columns
        auto& transactions = crel->transactions;
        transactions.erase(transactions.begin(), 
                upper_bound(transactions.begin(), transactions.end(), fstruct, TRSLessThan));
    }
    // then delete the transactions from the transaction history
    //gTransactionHistory.erase(gTransactionHistory.begin(), gTransactionHistory.upper_bound(f.transactionId));
#ifdef LPDEBUG
    LPTimer.forgets += LPTimer.getChrono(start);
#endif
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/////////////////// MAIN-READING STRUCTURES ///////////////////////
void ReaderTask(SRSWQueue<ReceivedMessage>& msgQ) {
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


int main(int argc, char**argv) {
    uint64_t numOfThreads = 1;
    if (argc > 1) {
        numOfThreads = strtol(argv[1], NULL, 10);
        if (numOfThreads == LONG_MAX)
            numOfThreads = 1;
    }
    cerr << "Number of threads: " << numOfThreads << endl;

    uint64_t MessageQSize  =50;
    uint64_t PendingMessages = MessageQSize-5;
    SRSWQueue<ReceivedMessage> msgQ(MessageQSize);
    std::thread readerTask(ReaderTask, std::ref(msgQ));
    
    SingleTaskPool workerThreads(numOfThreads, processPendingValidationsTask);
    workerThreads.initThreads();

try {
    //uint64_t msgs = 0;
    while (true) {
       
        // make sure the producer-ReaderTask- can write a new message
        if (gPendingTransactions.size() >= PendingMessages) {
            cerr << "it is all reserved: " << endl;
            // check if we have pending transactions to be processed
            processPendingMessages(workerThreads, msgQ);
            // release the space in msgQ
            cerr << "processed messagess" << endl;
        }
        
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
                Globals.state = GlobalState::VALIDATION;
                processValidationQueries(*reinterpret_cast<const ValidationQueries*>(msg.data.data())); 
                msgQ.registerDeq(res.refId);
                break;
            case MessageHead::Transaction: 
                Globals.state = GlobalState::TRANSACTION;
                //processTransaction(*reinterpret_cast<const Transaction*>(msg.data.data())); 
                //msgQ.registerDeq(res.refId);
                gPendingTransactions.push_back(res);
                // TODO - IMPORTANT DO NOT REGISTER THE DEQUEUE
                break;
            case MessageHead::Flush: { 
                // check if we have pending transactions to be processed
                processPendingMessages(workerThreads, msgQ);
                checkPendingValidations(workerThreads);
                Globals.state = GlobalState::FLUSH;
                processFlush(*reinterpret_cast<const Flush*>(msg.data.data())); 
                msgQ.registerDeq(res.refId);
                break;
                                     }
            case MessageHead::Forget: {
                // check if we have pending transactions to be processed
                processPendingMessages(workerThreads, msgQ);
                checkPendingValidations(workerThreads);
                Globals.state = GlobalState::FORGET;
                processForget(*reinterpret_cast<const Forget*>(msg.data.data())); 
                msgQ.registerDeq(res.refId);
                break;
                                      }
            case MessageHead::DefineSchema: 
                Globals.state = GlobalState::SCHEMA;
                processDefineSchema(*reinterpret_cast<const DefineSchema*>(msg.data.data()));
                gRelTransMutex = new std::mutex[gSchema.size()]();
                msgQ.registerDeq(res.refId);
                break;
            case MessageHead::Done: 
                {
#ifdef LPDEBUG
                    cerr << "  :::: " << LPTimer << endl; 
#endif              
                    delete[] gRelTransMutex;
                    readerTask.join();
                    workerThreads.destroy();
                    return 0;
                }
            default: cerr << "malformed message" << endl; abort(); // crude error handling, should never happen
        }

    }
} catch (const std::exception& e) { cerr <<  "exception " <<  e.what() << endl; }
}
//---------------------------------------------------------------------------


static uint64_t processPendingMessages(SingleTaskPool& pool, SRSWQueue<ReceivedMessage>& msgQ) {
    uint64_t msgs = 0;
    if (!gPendingTransactions.empty()) {
        msgs += processPendingTransactions(pool);
        for (auto& res : gPendingTransactions) msgQ.registerDeq(res.refId);
        gPendingTransactions.clear();
    }
    return msgs;
}

static std::atomic<uint64_t> gNextTrans;

static uint64_t processPendingTransactions(SingleTaskPool& pool) {
    if (gPendingTransactions.empty()) return 0;

    cerr << "ptrz: " << gPendingTransactions.size() << endl; 

    uint64_t ptrsz = gPendingTransactions.size();
    gNextTrans.store(0);
    pool.startAll(processPendingTransactionsTask);
    pool.waitAll();
    //gPendingTransactions.clear();

    return ptrsz;
}
static void processPendingTransactionsTask(uint32_t nThreads, uint32_t tid) {
    (void)nThreads; (void)tid;
    uint64_t ptrsz = gPendingTransactions.size();
    
    for (;true;) {
        uint64_t nextPos=gNextTrans++; 
        if (nextPos>=ptrsz) return;
        cerr << nextPos << ":" << ptrsz << endl;
        auto& res = gPendingTransactions[nextPos];
        
        processSingleTransaction(*reinterpret_cast<const Transaction*>(res.value->data.data()));
        //processSingleTransaction(*reinterpret_cast<const Transaction*>(gPendingTransactions[nextPos]->data.data()));
    }
}

static void processSingleTransaction(const Transaction& t) {
#ifdef LPDEBUG
    auto start = LPTimer.getChrono();
#endif
    const char* reader=t.operations;

    //cerr << t.transactionId << ":" << t.deleteCount << ":" << t.insertCount << endl;

    // Delete all indicated tuples
    for (uint32_t index=0;index!=t.deleteCount;++index) {
        auto& o=*reinterpret_cast<const TransactionOperationDelete*>(reader);
        auto& relation = gRelations[o.relationId];
        // TODO - lock here to make it to make all the deletions parallel naive locking first - 
        // TODO try to lock with try_lock and try again at the end if some relations failed
        {// start of lock_guard
        std::lock_guard<std::mutex> lk(gRelTransMutex[o.relationId]);
        for (const uint64_t* key=o.keys,*keyLimit=key+o.rowCount;key!=keyLimit;++key) {
            auto& rows = relation.insertedRows;
            auto lb = relation.insertedRows.find(*key);
            if (lb != rows.end()) {
                // update the relation transactions
                relation.transactions.push_back(move(TransStruct(t.transactionId, move(lb->second))));
                // remove the row from the relations table
                rows.erase(lb);
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
        uint64_t relCols = gSchema[o.relationId];
        // TODO - lock here to make it to make all the deletions parallel naive locking first - 
        // TODO try to lock with try_lock and try again at the end if some relations failed
        {// start of lock_guard
        std::lock_guard<std::mutex> lk(gRelTransMutex[o.relationId]);
        for (const uint64_t* values=o.values,*valuesLimit=values+(o.rowCount*relCols);values!=valuesLimit;values+=relCols) {
            vector<uint64_t> tuple(values, values+relCols);
            relation.transactions.push_back(move(TransStruct(t.transactionId, tuple)));
            relation.insertedRows[values[0]]=move(tuple);
        }
        }// end of lock_guard
        // advance to next Relation insertions
        reader+=sizeof(TransactionOperationInsert)+(sizeof(uint64_t)*o.rowCount*relCols);
    }

#ifdef LPDEBUG
    LPTimer.transactions += LPTimer.getChrono(start);
#endif
}




static uint64_t resIndexOffset = 0;
static std::atomic<uint64_t> gNextPending;

static void checkPendingValidations(SingleTaskPool &pool) {
    if (gPendingValidations.empty()) return;
    //if (!gPendingValidations.empty()) processPendingValidations(); 
#ifdef LPDEBUG
    auto start = LPTimer.getChrono();
#endif

    //cerr << gPendingValidations.size() << " " << endl;
    // find the min & max validation id
    // assuming that gPendingValidations is sorted on the validation Id
    resIndexOffset = gPendingValidations[0].validationId;   
    //auto gPVsz = gPendingValidations.size();
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
    //cerr << gPendingValidations.size() << " " << endl;
    //cerr << nThreads << ":" << tid << endl;

    uint64_t totalPending = gPendingValidations.size();
    //cerr << "start: " << rStart << " end: " << rEnd << endl;
    uint64_t resPos, vi;
    for (;true;) {

        // get a validation ID - atomic operation
        vi = gNextPending++;
        if (vi >= totalPending) { /*cerr << "exiting" << endl;*/  return; } // all pending finished
        //cerr << "got: " << vi << endl;

        auto& v = gPendingValidations[vi];
        resPos = v.validationId - resIndexOffset;
        auto& atoRes = gPendingResults[resPos];
        // check if someone else found a conflict already for this validation ID
        if (atoRes) continue;

        //cerr << "range: " << v.from << "-" << v.to << endl;
        TransStruct fromTRS(v.from);
        TransStruct toTRS(v.to);
        // TODO -  VERY NAIVE HERE - validate each query separately
        bool conflict = false, otherFinishedThis = false;
        for (auto& q : v.queries) {
            // avoid searching for the range of transactions too many times 
            auto& relation = gRelations[q.relationId];
            auto& transactions = relation.transactions;
            auto transFrom = std::lower_bound(transactions.begin(), transactions.end(), fromTRS, TRSLessThan);
            auto transTo = std::upper_bound(transactions.begin(), transactions.end(), toTRS, TRSLessThan);

            for(auto iter=transFrom; iter!=transTo; ++iter) {
                if (atoRes) { otherFinishedThis = true; /*cerr << "h" << endl;*/ break; }
                auto& tuple = iter->tuple;
                bool match=true;
                for (auto& c : q.predicates) {
                    // make the actual check
                    uint64_t tupleValue = tuple[c.column];
                    uint64_t queryValue=c.value;
                    bool result=false;
                    switch (c.op) {
                        case Query::Column::Equal: 
                            result=(tupleValue==queryValue); 
                            break;
                        case Query::Column::NotEqual: 
                            result=(tupleValue!=queryValue); 
                            break;
                        case Query::Column::Less: 
                            result=(tupleValue<queryValue); 
                            break;
                        case Query::Column::LessOrEqual: 
                            result=(tupleValue<=queryValue); 
                            break;
                        case Query::Column::Greater: 
                            result=(tupleValue>queryValue); 
                            break;
                        case Query::Column::GreaterOrEqual: 
                            result=(tupleValue>=queryValue); 
                            break;
                    } 
                    // there is one predicate not true so this whole query on this relation is false
                    if (!result) { match=false; break; }
                } // end of single query predicates
                if (match) { conflict=true; break; }    
            } // end of all transactions for this relation query
            if (conflict || otherFinishedThis) break;
        }
        // update the pending results to help other threads skip this validation 
        // if it has other parts
        if (conflict && !otherFinishedThis) atoRes = true;
    }
    /*
       {
       std::lock_guard<std::mutex> lk(gQueryResultsMutex);
       for (auto& p : localResults) gQueryResults[p.first]=p.second;
    //cerr << "global results: " << tid << endl;
    }
     */
}

