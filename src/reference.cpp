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
#include "BoundedSRSWQueue.hpp"
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
        columnCount = predicates.size();
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
bool QCSortOp (const Query::Column& left, const Query::Column& right) {
    if (left.op < right.op) return true;
    else if (right.op < left.op) return false;
    else if (left.column < right.column) return true;
    else if (right.column < left.column) return false;
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
bool LPQuerySizeLessThan(const LPQuery& left, const LPQuery& right) {
    if (left.columnCount < right.columnCount) return true;
    else if (right.columnCount < left.columnCount) return false;
    else return left.relationId < right.relationId;
}

struct LPValidation {
    uint64_t validationId;
    uint64_t from,to;
    vector<LPQuery> queries;

    LPValidation(const ValidationQueries& v, vector<LPQuery> q)
        : validationId(v.validationId), from(v.from), to(v.to), queries(move(q)) {}
};




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

static vector<uint32_t> gSchema;

static vector<LPValidation> gPendingValidations;

struct LPTimer_t {
    uint64_t validations;
    uint64_t validationsProcessing;
    uint64_t transactions;
    uint64_t flushes;
    uint64_t forgets;
    uint64_t reading;
    LPTimer_t() : validations(0), validationsProcessing(0), transactions(0), flushes(0), forgets(0), reading(0) {}
} LPTimer;
ostream& operator<< (ostream& os, const LPTimer_t& t) {
    os << "LPTimer [val: " << t.validations << " val-proc: " << t.validationsProcessing << " trans: " << t.transactions << " flush: " << t.flushes << " forget: " << t.forgets << " reads: " << t.reading << "]" << endl;
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
static void processPendingValidations();
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
    sort(queries.begin(), queries.end(), LPQuerySizeLessThan);

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

    // TODO - NEEDS A COMPLETE REFACTORING
    while ((!gQueryResults.empty())&&((*gQueryResults.begin()).first<=f.validationId)) {
        char c='0'+(*gQueryResults.begin()).second;
        cout.write(&c,1);
        gQueryResults.erase(gQueryResults.begin());
    }
    cout.flush();
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
// Read the message body and cast it to the desired type
template<typename Type> static const Type& readBody(istream& in,vector<char>& buffer,uint32_t len) {
#ifdef LPDEBUG
    auto start = getChrono();
#endif
    if (len > buffer.size()) buffer.resize(len);
    in.read(buffer.data(),len);
#ifdef LPDEBUG
    LPTimer.reading += getChrono(start);
#endif 
    return *reinterpret_cast<const Type*>(buffer.data());
}
//---------------------------------------------------------------------------

/////////////////// MAIN-READING STRUCTURES ///////////////////////

struct ReceivedMessage {
    MessageHead head;
    vector<char> data;
};
void ReaderTask(BoundedSRSWQueue<ReceivedMessage>& msgQ) {
    while (true) {
        // request place from the message queue - it blocks if full
        ReceivedMessage& msg = msgQ.reqNextEnq();
        auto& head = msg.head;
        auto& buffer = msg.data;
        // read the head of the message - type and len
        cin.read(reinterpret_cast<char*>(&head),sizeof(head));
        if (!cin) { cerr << "read error" << endl; abort(); } // crude error handling, should never happen

        if (head.type == MessageHead::Done) {
            msgQ.registerEnq(); 
            // exit the loop since the reader has finished its job
            return;        
        }

        // read the actual message content
        if (head.messageLen > buffer.size()) buffer.resize(head.messageLen);
        cin.read(buffer.data(), head.messageLen);
        msgQ.registerEnq();
    }
    return;
}

class Barrier
{
    private:
        std::mutex _mutex;
        std::condition_variable _cv;
        std::size_t _count;
        uint32_t mSize;
    public:
        explicit Barrier(std::size_t count) : _count{count}, mSize(count) { }
        void wait()
        {
            std::unique_lock<std::mutex> lock{_mutex};
            if (--_count == 0) {
                _cv.notify_all();
            } else {
                _cv.wait(lock, [this] { return _count == 0; });
            }
        }

        void reset() {
            std::lock_guard<std::mutex> lock(_mutex);
            _count = mSize;
        }
};

class SingleTaskPool {
    public:
        SingleTaskPool(uint32_t tsz) : mNumOfThreads(tsz), mPoolStopped(false), 
            mActiveThreadsReq(0), mBarrier(tsz+1) {
        }

        template<typename Func>
        void initThreads(Func func) {
            for (uint32_t i=0; i<mNumOfThreads; ++i) {
                mThreads.push_back(
                        std::thread([i,func,this](){
                            uint32_t tid = i;
                            //cerr << "tid:" << tid << endl;
                            while(true) {
                            std::unique_lock<std::mutex> lk(mMutex);
                            mCondActive.wait(lk, [this]{return mActiveThreadsReq > 0 || mPoolStopped;});
                            --mActiveThreadsReq;
                            lk.unlock();

                            // check if the pool is still active
                            if (mPoolStopped) return;
                            // run the function passing the thread ID
                            func(tid);
                            // wait after the execution for the other threads
                            mBarrier.wait();
                            }
                            }));
            }
        }

        void startAll() {
            std::lock_guard<std::mutex> lk(mMutex);
            mActiveThreadsReq = mNumOfThreads;
            mBarrier.reset();
            mCondActive.notify_all();
        }

        void waitAll() {
            mBarrier.wait();
        }

        void destroy() {
            std::lock_guard<std::mutex> lk(mMutex);
            mPoolStopped = true;
            mActiveThreadsReq = 0;
            mCondActive.notify_all();
            for(auto& t : mThreads) t.join();
        }

    private:
        uint32_t mNumOfThreads;
        vector<std::thread> mThreads;
        bool mPoolStopped;
        uint32_t mActiveThreadsReq;
        std::condition_variable mCondActive;
        std::mutex mMutex;

        Barrier mBarrier;
};

void dummy(uint32_t tid) {
    cerr << "thread executing : " << tid << endl;
}

int main(int argc, char**argv) {
    // TODO - MAYBE some optimizations on reading
    ios::sync_with_stdio(false);

    uint64_t numOfThreads = 1;
    if (argc > 1) {
        numOfThreads = strtol(argv[1], NULL, 10);
        if (numOfThreads == LONG_MAX)
            numOfThreads = 1;
    }

    cerr << "Number of threads: " << numOfThreads << endl;

    BoundedSRSWQueue<ReceivedMessage> msgQ(10000);

    SingleTaskPool validationThreads(numOfThreads);
    validationThreads.initThreads(dummy);

    std::thread readerTask(ReaderTask, std::ref(msgQ));

    //vector<char> message;
    //MessageHead head;
    while (true) {
        // Retrieve the message
        ReceivedMessage& msg = msgQ.reqNextDeq();
        auto& head = msg.head;
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
                    cerr << "\ttimings:: " << LPTimer << endl; 
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

static inline void checkPendingValidations(SingleTaskPool &pool) {
    if (gPendingValidations.empty()) return;
#ifdef LPDEBUG
    auto start = getChrono();
#endif

    pool.startAll();
    processPendingValidations();
    pool.waitAll();

#ifdef LPDEBUG
    LPTimer.validationsProcessing += getChrono(start);
#endif
}

static void processPendingValidations() {
    //cerr << gPendingValidations.size() << endl;
    uint64_t vsz=gPendingValidations.size();
    for (uint64_t vi=0; vi<vsz; ++vi) {
        auto& v = gPendingValidations[vi];
        //cerr << "range: " << v.from << "-" << v.to << endl;

        TransStruct fromTRS(v.from);
        TransStruct toTRS(v.to);
        // TODO -  VERY NAIVE HERE - validate each query separately
        bool conflict=false;
        for (auto& q : v.queries) {
            //auto& q=v.queries[index];
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
        gQueryResults[v.validationId]=conflict;
    }

    gPendingValidations.clear();

}
