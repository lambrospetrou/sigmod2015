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

typedef vector<uint64_t> tuple_t;

// Custom data structures to hold data
struct CTransStruct {
    uint64_t trans_id;
    uint64_t value;
    tuple_t *tuple;
    CTransStruct (uint64_t tid, uint64_t v, tuple_t *t) : trans_id(tid), value(v), tuple(t) {}
    bool operator< (const CTransStruct& o) { 
        if (value < o.value) return true;
        else if (o.value < value) return false;
        else return trans_id < o.trans_id;
    }
    static bool CompValTrans(const CTransStruct& l, const CTransStruct& r) {
        if (l.value < r.value) return true;
        else if (r.value < l.value) return false;
        else return l.trans_id < r.trans_id;
    }
    static bool CompTransVal(const CTransStruct& l, const CTransStruct& r) {
        if (l.trans_id < r.trans_id) return true;
        else if (r.trans_id < l.trans_id) return false;
        else return l.value < r.value;
    }
    static bool CompValOnly(const CTransStruct& l, const CTransStruct& r) {
        return l.value < r.value;
    }
};
std::ostream& operator<< (std::ostream& os, const CTransStruct& o) {
    os << "{" << o.trans_id << "-" << o.value << "}";
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
    //std::map<uint64_t, vector<CTransStruct>> transactions;
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
struct RelationStruct {
    vector<TransStruct> transactions;
    unordered_map<uint32_t, std::unique_ptr<tuple_t>> insertedRows;
};
static vector<RelationStruct> gRelations;

struct RelationColumns {
    vector<ColumnStruct> columns;
};
static vector<RelationColumns> gRelColumns;


struct TransOperation {
    uint32_t rel_id;
    tuple_t *tuple;
    TransOperation(uint32_t relid, tuple_t *t):
        rel_id(relid), tuple(t) {}
};
struct TransactionStruct {
    uint64_t trans_id;
    vector<TransOperation> operations;
    TransactionStruct(uint64_t tid, vector<TransOperation> ops) : trans_id(tid), operations(move(ops)) {}  
};
static vector<TransactionStruct> gTransactionHistory;

static vector<uint32_t> gSchema;

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


///--------------------------------------------------------------------------
///--------------------------------------------------------------------------

static void processDefineSchema(const DefineSchema& d) {
    gSchema.clear();
    gSchema.insert(gSchema.begin(),d.columnCounts,d.columnCounts+d.relationCount);

    gRelations.clear();
    gRelations.resize(d.relationCount);
    gRelColumns.resize(d.relationCount);
    for(uint32_t ci=0; ci<d.relationCount; ++ci) {
        //gRelations[ci].columns.resize(gSchema[ci]);
        gRelColumns[ci].columns.resize(gSchema[ci]);
    }
}
//---------------------------------------------------------------------------

static void processTransaction(const Transaction& t, const vector<char>& tdata) {
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
    //static uint64_t bSize = 5000;
    {
        std::lock_guard<std::mutex> lk(gPendingValidationsMutex);
        /*
           while (trRange > bSize) {
        //cerr << batchPos << "-" << batchPos+bSize-1 << endl;
        gPendingValidations.push_back(move(LPValidation(v.validationId, batchPos, batchPos+bSize-1, queries)));    
        trRange -= bSize; batchPos += bSize;
        }
         */
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
    (void)f.transactionId;
    /*
    // Delete the transactions from inside the columns in the relations
    for(auto crel=gRelations.begin(), iend=gRelations.end(); crel!=iend; ++crel) {
    // delete this transaction from the lastRel columns
    auto& transactions = crel->transactions;
    transactions.erase(transactions.begin(), 
    upper_bound(transactions.begin(), transactions.end(), f.transactionId, TRSLessThan));
    }
     */
    // then delete the transactions from the transaction history
    /*
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

void parseValidation(uint32_t nThreads, uint32_t tid, void *args) {
    (void)tid; (void)nThreads;
    ParseValidationStruct *pvs = static_cast<ParseValidationStruct*>(args);
    processValidationQueries(*reinterpret_cast<const ValidationQueries*>(pvs->msg->data.data()), pvs->msg->data); 
    pvs->msgQ->registerDeq(pvs->refId);
    pvs->memQ->free(pvs->memRefId);
}

static uint64_t gTotalTransactions = 0, gTotalTuples = 0;

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

    MultiTaskPool multiPool(numOfThreads-2);
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
                        ++gTotalValidations; // this is just to count the total validations....not really needed!

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
                    Globals.state = GlobalState::TRANSACTION;
                    processTransaction(*reinterpret_cast<const Transaction*>(msg.data.data()), msg.data); 
                    msgQ.registerDeq(res.refId);
                    break;
                case MessageHead::Flush:  
                    // check if we have pending transactions to be processed
                    multiPool.waitAll();
                    checkPendingValidations(workerThreads);
                    Globals.state = GlobalState::FLUSH;
                    processFlush(*reinterpret_cast<const Flush*>(msg.data.data())); 
                    msgQ.registerDeq(res.refId);
                    break;

                case MessageHead::Forget: 
                    // check if we have pending transactions to be processed
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
#endif
    const char* reader=t.operations;
    ++gTotalTransactions; 
    //cerr << t.transactionId << ":" << t.deleteCount << ":" << t.insertCount << endl;

    //gTransactionHistory.push_back(move(TransactionStruct(t.transactionId, vector<TransOperation>())));
    //vector<TransOperation>& operations = gTransactionHistory.back().operations;

    // Delete all indicated tuples
    for (uint32_t index=0;index!=t.deleteCount;++index) {
        auto& o=*reinterpret_cast<const TransactionOperationDelete*>(reader);
        auto& relation = gRelations[o.relationId];
        auto& relColumns = gRelColumns[o.relationId].columns;
        uint32_t relCols = gSchema[o.relationId];

        //cerr << " " << o.relationId;

        // TODO - lock here to make it to make all the deletions parallel naive locking first - 
        // TODO try to lock with try_lock and try again at the end if some relations failed
        {// start of lock_guard
            //std::lock_guard<std::mutex> lk(gRelTransMutex[o.relationId]);


            for (uint32_t col=0; col<relCols; ++col) {
                relColumns[col].transactions.push_back(move(std::make_pair(t.transactionId, move(vector<CTransStruct>()))));
            }

            for (const uint64_t* key=o.keys,*keyLimit=key+o.rowCount;key!=keyLimit;++key) {
                auto& rows = relation.insertedRows;
                auto lb = relation.insertedRows.find(*key);
                if (lb != rows.end()) {
                    // update the relation transactions - transfer ownership of the tuple
                    relation.transactions.push_back(move(TransStruct(t.transactionId, move(lb->second))));

                    // insert into transaction history
                    tuple_t *tpl = relation.transactions.back().tuple.get();
                    //operations.push_back(move(TransOperation(o.relationId, tpl)));
                    
//                   cerr << *tpl << endl;


                    for (uint32_t col=0; col<relCols; ++col) {
                        relColumns[col].transactions.back().second.push_back(move(CTransStruct(t.transactionId, (*tpl)[col], tpl)));
                    }

                    // insert the tuples only in column 0 - primary key
                    //relation.columns[0].transactions[(*tpl)[0]].push_back(move(CTransStruct(t.transactionId, (*tpl)[0], tpl)));
                    //gRelColumns[o.relationId].columns[0].transactions.push_back(move(CTransStruct(t.transactionId, (*tpl)[0], tpl)));

                    // remove the row from the relations table
                    rows.erase(lb);
                    ++gTotalTuples;
                }
            }

            // SORT THE VALUES
            for (uint32_t col=0; col<relCols; ++col) {
                auto& vec = relColumns[col].transactions.back().second;
                std::sort(vec.begin(), vec.end(), CTransStruct::CompValOnly);
            }
/*
            for (auto& ctr : relColumns[2].transactions.back().second) {
                cerr << *ctr.tuple << endl;
            }
            cerr << endl << endl;
*/
        }// end of lock_guard

        // advance to the next Relation deletions
        reader+=sizeof(TransactionOperationDelete)+(sizeof(uint64_t)*o.rowCount);
    }

    //cerr << endl;

    // Insert new tuples
    for (uint32_t index=0;index!=t.insertCount;++index) {
        auto& o=*reinterpret_cast<const TransactionOperationInsert*>(reader);
        auto& relation = gRelations[o.relationId];
        auto& relColumns = gRelColumns[o.relationId].columns;
        uint32_t relCols = gSchema[o.relationId];
        // TODO - lock here to make it to make all the deletions parallel naive locking first - 
        // TODO try to lock with try_lock and try again at the end if some relations failed
        {// start of lock_guard
            //std::lock_guard<std::mutex> lk(gRelTransMutex[o.relationId]);

            // TODO change this in order to use the same vectors for DELETES & INSERTS
            for (uint32_t col=0; col<relCols; ++col) {
                relColumns[col].transactions.push_back(move(std::make_pair(t.transactionId, move(vector<CTransStruct>()))));
            }

            for (const uint64_t* values=o.values,*valuesLimit=values+(o.rowCount*relCols);values!=valuesLimit;values+=relCols) {
                unique_ptr<tuple_t> tptr(new tuple_t(values, values+relCols));
                relation.transactions.push_back(move(TransStruct(t.transactionId, move(tptr))));
                unique_ptr<tuple_t> tptr2(new tuple_t(values, values+relCols));
                relation.insertedRows[values[0]]=move(tptr2);

                // history holds RAW pointers to the tuple
                tuple_t *tpl = relation.transactions.back().tuple.get();
                //operations.push_back(move(TransOperation(o.relationId, tpl)));

                // insert the tuples only in column 0 - primary key
                //relation.columns[0].transactions[values[0]].push_back(move(CTransStruct(t.transactionId, values[0], tpl)));
                //gRelColumns[o.relationId].columns[0].transactions.push_back(move(CTransStruct(t.transactionId, (*tpl)[0], tpl)));
                for (uint32_t col=0; col<relCols; ++col) {
                    relColumns[col].transactions.back().second.push_back(move(CTransStruct(t.transactionId, values[col], tpl)));
                }

                ++gTotalTuples;
            }

            // SORT THE VALUES
            for (uint32_t col=0; col<relCols; ++col) {
                auto& vec = relColumns[col].transactions.back().second;
                std::sort(vec.begin(), vec.end(), CTransStruct::CompValOnly);
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
    //resIndexOffset = gPendingValidations[0].validationId;   
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

        //if (v.validationId == 19631)  cerr << v.queries << endl;
        //if (v.validationId == 6673)  cerr << v.queries << endl;
        //if (v.validationId == 4)  cerr << v.queries << endl;

        // TODO -  VERY NAIVE HERE - validate each query separately
        bool conflict = false, otherFinishedThis = false;
        for (auto& q : v.queries) {
            if (atoRes) { otherFinishedThis = true; /*cerr << "h" << endl;*/ break; }

            auto& relColumns = gRelColumns[q.relationId].columns;
            
            // avoid searching for the range of transactions too many times 
            //auto& relation = gRelations[q.relationId];
            //auto& transactions = relation.transactions;
            //auto transFrom = std::lower_bound(transactions.begin(), transactions.end(), v.from, TRSLessThan);
            //auto transTo = std::upper_bound(transactions.begin(), transactions.end(), v.to, TRSLessThan);
            
            //cerr << ":: validation: " << v.validationId << " query: " << q << endl;
            // protect from the case where there is no single predicate
            //if (q.predicates.empty()) { cerr << "empty: " << v.validationId << endl; conflict = true; break; }
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
            //auto transFrom = transactions.begin();
            //auto transTo = transactions.end();

            //cerr << "before: " << (transTo-transFrom) << " for col: " << pFirst.column << "-" << pFirst.value << endl;

            // take only the transactions that we are interested in
            /*
            while (transFrom != transactions.end() && transFrom->first < v.from) ++transFrom;
            transTo = transFrom;
            while (transTo != transactions.end() && transTo->first <= v.to) ++transTo;
            */

            //if (v.validationId == 4)
            //cerr << "after: " << v.from << "-" << v.to << "=" << (transTo-transFrom) << " for col: " << pFirst.column << "-" << pFirst.value << endl;

            for(auto iter=transFrom; iter!=transTo; ++iter) {  
                // TODO - check for transactions
                //if (transFrom->first < v.from || transFrom->first > v.to) { cerr << "wrong trans" << endl; continue;}

                auto& transValues = iter->second;
                uint32_t pFrom = 0;
                decltype(transValues.begin()) tupFrom, tupTo;
                // find the valid tuples using range binary searches based on the first predicate
                if (pFirst.op == Query::Column::Equal) {
                    tupFrom = std::lower_bound(transValues.begin(), transValues.end(), pFirst.value, CTRSValueLessThan);                    
                    tupTo = std::upper_bound(transValues.begin(), transValues.end(), pFirst.value, CTRSValueLessThan);                   
                    pFrom = 1;
                } else if (pFirst.op == Query::Column::Less) {
                    tupFrom = transValues.begin();                    
                    tupTo = std::lower_bound(transValues.begin(), transValues.end(), pFirst.value, CTRSValueLessThan);                   
                    pFrom = 1;
                } else if (pFirst.op == Query::Column::LessOrEqual) {
                    tupFrom = transValues.begin();                    
                    tupTo = std::upper_bound(transValues.begin(), transValues.end(), pFirst.value, CTRSValueLessThan);                   
                    pFrom = 1;
                } else if (pFirst.op == Query::Column::Greater) {
                    tupFrom = std::upper_bound(transValues.begin(), transValues.end(), pFirst.value, CTRSValueLessThan);                    
                    tupTo = transValues.end();                   
                    pFrom = 1;
                } else if (pFirst.op == Query::Column::GreaterOrEqual) {
                    tupFrom = std::lower_bound(transValues.begin(), transValues.end(), pFirst.value, CTRSValueLessThan);                    
                    tupTo = transValues.end();                   
                    pFrom = 1;
                } else {
                    tupFrom = transValues.begin();
                    tupTo = transValues.end();
                }
                
                //cerr << "tup diff " << (tupTo - tupFrom) << endl; 

                //cerr << "first tuple " << *tupFrom->tuple << endl;

                //if (v.validationId == 4)
                //cerr << "trans " << iter->first << endl;

                for(; tupFrom!=tupTo; ++tupFrom) {  
                    tuple_t& tuple = *tupFrom->tuple;
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

