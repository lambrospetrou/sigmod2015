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
#include <memory>
#include "BoundedSRSWQueue.hpp"
#include "LockFreeBoundedSRSWQueue.hpp"
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
    LPQuery(const Query& q) : relationId(q.relationId), columnCount(q.columnCount), 
    predicates(vector<Query::Column>(q.columnCount)) {
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

static map<uint64_t, TransactionStruct> gTransactionHistory;
static map<uint64_t,bool> gQueryResults;
static std::mutex gQueryResultsMutex;

static vector<uint32_t> gSchema;


template <typename T>
class atomic_wrapper {
    private:
        std::atomic<T> _a;
    public:
        atomic_wrapper()
            :_a() {}

        atomic_wrapper(const std::atomic<T> &a)
            :_a(a.load()){}

        atomic_wrapper(const atomic_wrapper &other)
            :_a(other._a.load()){}

        atomic_wrapper &operator=(const atomic_wrapper &other) {
            _a.store(other._a.load());
        }

        void store(T desired) {_a.store(desired);}
        void store(T desired) volatile {_a.store(desired);}
        T load() const { return _a.load(); }
        T load() const volatile { return _a.load(); }
};
struct LPValidation {
    uint64_t validationId;
    uint64_t from,to;
    vector<LPQuery> queries;
    LPValidation(const ValidationQueries& v, vector<LPQuery> q)
        : validationId(v.validationId), from(v.from), to(v.to), queries(move(q)) {}
};
static vector<LPValidation> gPendingValidations;
//static unique_ptr<std::atomic<bool>[]> gPendingResults;
static vector<atomic_wrapper<bool>> gPendingResults;

struct LPTimer_t {
    uint64_t validations;
    uint64_t validationsProcessing;
    uint64_t transactions;
    uint64_t flushes;
    uint64_t forgets;
    uint64_t reading;
    uint64_t readingTotal;
    LPTimer_t() : validations(0), validationsProcessing(0), transactions(0), flushes(0), forgets(0), reading(0), readingTotal(0) {}
} LPTimer;
ostream& operator<< (ostream& os, const LPTimer_t& t) {
    os << "LPTimer [val: " << t.validations << " val-proc: " << t.validationsProcessing << " trans: " << t.transactions << " flush: " << t.flushes << " forget: " << t.forgets << " reads: " << t.reading << "/" << t.readingTotal<< "]" << endl;
    return os;
}
uint64_t getChronoMicro() {      
    // return nanoseconds
    struct timeval start;
    gettimeofday(&start, NULL);
    // tv_sec = seconds | tv_usecs = microseconds
    return (start.tv_sec * 1000000LL) + start.tv_usec;
}
uint64_t getChronoMicro(uint64_t start) {      
    return getChronoMicro() - start;
}
uint64_t getChrono() {
    return getChronoMicro();
}
uint64_t getChrono(uint64_t start) {
    return getChronoMicro() - start;
}

struct GlobalState {
    enum State : uint32_t { SCHEMA, TRANSACTION, VALIDATION, FORGET, FLUSH };
    State state;

} Globals;

//---------------------------------------------------------------------------

// JUST SOME FUNCTION DECLARATIONS THAT ARE DEFINED BELOW
class SingleTaskPool;
static inline void checkPendingValidations(SingleTaskPool&);

template<class T> using SRSWQueue = LockFreeBoundedSRSWQueue<T>;


///--------------------------------------------------------------------------
///--------------------------------------------------------------------------

static void processDefineSchema(const DefineSchema& d) {
    gSchema.clear();
    gSchema.insert(gSchema.begin(),d.columnCounts,d.columnCounts+d.relationCount);

    gRelations.clear();
    // TODO - maybe always allocate the maximum of 1000 columns to avoid the following resize loop
    gRelations.resize(d.relationCount);
    for(uint32_t ci=0; ci<d.relationCount; ++ci) {
        gRelations[ci].columns.resize(gSchema[ci]);
    }
}
//---------------------------------------------------------------------------
static void processTransaction(const Transaction& t) {
#ifdef LPDEBUG
    auto start = getChrono();
#endif
    //vector<TransOperation> operations;
    const char* reader=t.operations;

    // Delete all indicated tuples
    for (uint32_t index=0;index!=t.deleteCount;++index) {
        auto& o=*reinterpret_cast<const TransactionOperationDelete*>(reader);
        for (const uint64_t* key=o.keys,*keyLimit=key+o.rowCount;key!=keyLimit;++key) {
            auto& relation = gRelations[o.relationId];
            auto& rows = relation.insertedRows;
            auto lb = relation.insertedRows.find(*key);
            if (lb != rows.end()) {
                //uint64_t rowid = lb->second[0];
                // update the relation transactions
                relation.transactions.push_back(move(TransStruct(t.transactionId, move(lb->second))));
                // update the transactions history record - so far we to not need tuple
                //operations.push_back(TransOperation(o.relationId, rowid, relation.transactions.back().tuple));
                //operations.push_back(move(TransOperation(o.relationId, rowid)));
                // remove the row from the relations table
                rows.erase(lb);
            }
        }
        reader+=sizeof(TransactionOperationDelete)+(sizeof(uint64_t)*o.rowCount);
    }

    // Insert new tuples
    for (uint32_t index=0;index!=t.insertCount;++index) {
        auto& o=*reinterpret_cast<const TransactionOperationInsert*>(reader);
        for (const uint64_t* values=o.values,*valuesLimit=values+(o.rowCount*gSchema[o.relationId]);values!=valuesLimit;values+=gSchema[o.relationId]) {
            vector<uint64_t> tuple(gSchema[o.relationId]);
            //tuple.insert(tuple.begin(),values,values+gSchema[o.relationId]);
            memcpy(tuple.data(),values,gSchema[o.relationId]*sizeof(uint64_t));

            // add the tuple to this transaction operations and to the relations table
            //operations.push_back(move(TransOperation(o.relationId, tuple[0])));

            gRelations[o.relationId].transactions.push_back(move(TransStruct(t.transactionId, tuple)));
            gRelations[o.relationId].insertedRows[values[0]]=move(tuple);
        }
        reader+=sizeof(TransactionOperationInsert)+(sizeof(uint64_t)*o.rowCount*gSchema[o.relationId]);
    }

    //gTransactionHistory.insert(move(std::pair<uint64_t, TransactionStruct>(t.transactionId, TransactionStruct(move(operations))))); 

#ifdef LPDEBUG
    LPTimer.transactions += getChrono(start);
#endif
}
//---------------------------------------------------------------------------
static void processValidationQueries(const ValidationQueries& v) {
#ifdef LPDEBUG
    auto start = getChrono();    
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

    gPendingValidations.push_back(move(LPValidation(v, move(queries))));

    //cerr << "Success Validate: " << v.validationId << " ::" << v.from << ":" << v.to << " result: " << conflict << endl;
#ifdef LPDEBUG
    LPTimer.validations += getChrono(start);
#endif
}
//---------------------------------------------------------------------------
static void processFlush(const Flush& f) {
#ifdef LPDEBUG
    auto start = getChrono();
#endif
    char zero = 48;
    char one = 49;
    // TODO - NEEDS A COMPLETE REFACTORING
    while ((!gQueryResults.empty())&&((*gQueryResults.begin()).first<=f.validationId)) {
        cout << (gQueryResults.begin()->second == true ? one : zero); 
        //cerr << (gQueryResults.begin()->second == true ? one : zero) << " "; 
        gQueryResults.erase(gQueryResults.begin());
    }
    cout.flush();

    //cerr << "finished flushing" << endl;
#ifdef LPDEBUG
    LPTimer.flushes += getChrono(start);
#endif
}
//---------------------------------------------------------------------------

static void processForget(const Forget& f) {
#ifdef LPDEBUG
    //cerr << "Forget: " << f.transactionId << endl;
    auto start = getChrono();
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
    gTransactionHistory.erase(gTransactionHistory.begin(), gTransactionHistory.upper_bound(f.transactionId));
#ifdef LPDEBUG
    LPTimer.forgets += getChrono(start);
#endif
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/////////////////// MAIN-READING STRUCTURES ///////////////////////
struct ReceivedMessage {
    MessageHead head;
    vector<char> data;
};
void ReaderTask(SRSWQueue<ReceivedMessage>& msgQ) {
#ifdef LPDEBUG
    auto start = getChrono();
#endif
    while (true) {
        // request place from the message queue - it blocks if full
        ReceivedMessage& msg = msgQ.reqNextEnq();
        auto& head = msg.head;
        auto& buffer = msg.data;
        // read the head of the message - type and len
        // Read the message body and cast it to the desired type
        cin.read(reinterpret_cast<char*>(&head),sizeof(head));
#ifdef LPDEBUG // I put the inner timer here to avoid stalls in the msgQ
        auto startInner = getChrono();
#endif
        if (!cin) { cerr << "read error" << endl; abort(); } // crude error handling, should never happen

        if (head.type == MessageHead::Done) {
            msgQ.registerEnq(); 
            // exit the loop since the reader has finished its job
            break;        
        }

        // read the actual message content
        if (head.messageLen > buffer.size()) buffer.resize(head.messageLen);
        cin.read(buffer.data(), head.messageLen);
#ifdef LPDEBUG
        LPTimer.reading += getChrono(startInner);
#endif 
        msgQ.registerEnq();
    }
#ifdef LPDEBUG
    LPTimer.readingTotal += getChrono(start);
#endif 
    return;
}

static void processPendingValidationsTask(uint32_t nThreads, uint32_t tid);

int main(int argc, char**argv) {
    uint64_t numOfThreads = 1;
    if (argc > 1) {
        numOfThreads = strtol(argv[1], NULL, 10);
        if (numOfThreads == LONG_MAX)
            numOfThreads = 1;
    }

    cerr << "Number of threads: " << numOfThreads << endl;

    SingleTaskPool validationThreads(numOfThreads, processPendingValidationsTask);
    validationThreads.initThreads();

    SRSWQueue<ReceivedMessage> msgQ(100);
    std::thread readerTask(ReaderTask, std::ref(msgQ));

    while (true) {
        // Retrieve the message
        //cerr << "try for incoming" << endl;
        ReceivedMessage& msg = msgQ.reqNextDeq();
        auto& head = msg.head;
        //cerr << "incoming: " << head.type << endl;
        // And interpret it
        switch (head.type) {
            case MessageHead::ValidationQueries: 
                Globals.state = GlobalState::VALIDATION;
                processValidationQueries(*reinterpret_cast<const ValidationQueries*>(msg.data.data())); 
                break;
            case MessageHead::Transaction: 
                Globals.state = GlobalState::TRANSACTION;
                processTransaction(*reinterpret_cast<const Transaction*>(msg.data.data())); 
                break;
            case MessageHead::Flush: 
                checkPendingValidations(validationThreads);
                Globals.state = GlobalState::FLUSH;
                processFlush(*reinterpret_cast<const Flush*>(msg.data.data())); 
                break;
            case MessageHead::Forget: 
                checkPendingValidations(validationThreads);
                Globals.state = GlobalState::FORGET;
                processForget(*reinterpret_cast<const Forget*>(msg.data.data())); 
                break;
            case MessageHead::DefineSchema: 
                Globals.state = GlobalState::SCHEMA;
                processDefineSchema(*reinterpret_cast<const DefineSchema*>(msg.data.data())); 
                break;
            case MessageHead::Done: 
                {
#ifdef LPDEBUG
                    cerr << "  :::: " << LPTimer << endl; 
#endif              
                    readerTask.join();
                    validationThreads.destroy();
                    return 0;
                }
            default: cerr << "malformed message" << endl; abort(); // crude error handling, should never happen
        }

        msgQ.registerDeq();
    }
}
//---------------------------------------------------------------------------

static uint64_t resIndexOffset = 0;
static std::atomic<uint64_t> gNextPending;

static inline void checkPendingValidations(SingleTaskPool &pool) {
    if (gPendingValidations.empty()) return;
    //if (!gPendingValidations.empty()) processPendingValidations(); 
#ifdef LPDEBUG
    auto start = getChrono();
#endif

    //cerr << gPendingValidations.size() << " " << endl;
    // find the min & max validation id
    // assuming that gPendingValidations is sorted on the validation Id
    resIndexOffset = gPendingValidations[0].validationId;    
    if (gPendingValidations.size() > gPendingResults.size())
        gPendingResults.resize(gPendingValidations.size());
    memset(gPendingResults.data(), 0, sizeof(atomic_wrapper<bool>)*gPendingResults.size());
    //gPendingResults.reset(new std::atomic<bool>[gPendingValidations.size()]());
    //    for (uint64_t i=0; i<gPendingValidations.size();++i) if (gPendingResults[i].load()) cerr << "ERRORORROROROOROROROROROR" << endl;
    gNextPending.store(0);

    //cerr << "before start!" << endl;
    pool.startAll();
    //processPendingValidations();
    //cerr << "before wait!" << endl;
    pool.waitAll();
    //cerr << "after wait!" << endl;

    // update the results - you can get the validation id by adding resIndexOffset to the position
    for (uint64_t i=0, sz=gPendingValidations.size(); i<sz; ++i) { 
        //cerr << i << ":" << gPendingResults[i].load() << " ";
        gQueryResults[gPendingValidations[i].validationId] = gPendingResults[i].load();
    }
    gPendingValidations.clear();

#ifdef LPDEBUG
    LPTimer.validationsProcessing += getChrono(start);
#endif
}

static void processPendingValidationsTask(uint32_t nThreads, uint32_t tid) {
    (void)tid; (void)nThreads;// to avoid unused warning
    //cerr << gPendingValidations.size() << " " << endl;
    //vector<pair<uint64_t, bool>> localResults;
    //cerr << nThreads << ":" << tid << endl;

    //uint64_t vsz=gPendingValidations.size();
    //uint64_t batchSize = vsz/nThreads;
    //uint64_t remSize = vsz%nThreads;
    //uint64_t rStart = tid*batchSize + (tid < remSize ? tid : remSize);
    // do not give all the remaining to the last node but each thread take sone more
    //uint64_t rEnd = rStart + batchSize + (tid < remSize);
    uint64_t totalPending = gPendingValidations.size();
    //cerr << "start: " << rStart << " end: " << rEnd << endl;
    //for (uint64_t vi=rStart; vi<rEnd; ++vi) {
    for (;true;) {

        // get a validation ID - atomic operation
        uint64_t vi = gNextPending++;
        if (vi >= totalPending) { /*cerr << "exiting" << endl;*/  return; } // all pending finished
        //cerr << "got: " << vi << endl;

        auto& v = gPendingValidations[vi];
        auto resPos = v.validationId - resIndexOffset;
        auto& atoRes = gPendingResults[resPos];
        // check if someone else found a conflict already for this validation ID
        if (atoRes.load()) continue;

        //cerr << "range: " << v.from << "-" << v.to << endl;
        TransStruct fromTRS(v.from);
        TransStruct toTRS(v.to);
        // TODO -  VERY NAIVE HERE - validate each query separately
        bool conflict = false, otherFinishedThis = false;
        for (auto& q : v.queries) {
            if (atoRes.load()) { otherFinishedThis = true; break; }
            // avoid searching for the range of transactions too many times 
            auto& relation = gRelations[q.relationId];
            auto& transactions = relation.transactions;
            auto transFrom = std::lower_bound(transactions.begin(), transactions.end(), fromTRS, TRSLessThan);
            auto transTo = std::upper_bound(transactions.begin(), transactions.end(), toTRS, TRSLessThan);

            for(auto iter=transFrom; iter!=transTo; ++iter) {
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
            if (conflict) break;
        }
        // update the pending results to help other threads skip this validation 
        // if it has other parts
        if (conflict && !otherFinishedThis) atoRes.store(true);
    }
    /*
       {
       std::lock_guard<std::mutex> lk(gQueryResultsMutex);
       for (auto& p : localResults) gQueryResults[p.first]=p.second;
    //cerr << "global results: " << tid << endl;
    }
     */
}

