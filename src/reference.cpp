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
#include <array>
#include <cassert>
#include <cstdint>
#include <algorithm>
#include <utility>
#include <chrono>
#include <cstring>
//---------------------------------------------------------------------------
using namespace std;
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
    vector<Query::Column> columns;
    LPQuery() : relationId(-1), columnCount(0), columns(vector<Query::Column>()) {}
    LPQuery(const Query& q) : relationId(q.relationId), columnCount(q.columnCount), 
            columns(vector<Query::Column>(q.columnCount)) {
        memcpy(columns.data(), q.columns, sizeof(Query::Column)*columnCount);
        std::stable_sort(columns.begin(), columns.end());
        columns.resize(std::distance(columns.begin(), std::unique(columns.begin(), columns.end())));
        //if (columns.size() != columnCount) cerr << "diff: " << columnCount-columns.size() << endl;
        columnCount = columns.size();
    }
};
bool operator< (const LPQuery& left, const LPQuery& right) {
    return left.relationId < right.relationId;
}
bool operator== (const LPQuery& left, const LPQuery& right)  {
    if (left.relationId != right.relationId) return false;
    if (left.columnCount != right.columnCount) return false;
    return left.columns == right.columns;
}
ostream& operator<< (ostream& os, const LPQuery& o) {
    os << "{" << o.relationId << "-" << o.columnCount << ":: " << o.columns << "}";
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
    vector<uint64_t> tuple;
    TransOperation(uint32_t relid, uint64_t rowid, vector<uint64_t> t):
        rel_id(relid), row_id(rowid), tuple(t) {}
    TransOperation(uint32_t relid, uint64_t rowid):
        rel_id(relid), row_id(rowid), tuple(vector<uint64_t>()) {}
};
struct TransactionStruct {
    //    uint64_t trans_id; // we will use direct indexing
    vector<TransOperation> operations;
    TransactionStruct(vector<TransOperation> ops) : operations(ops) {}
};
static map<uint64_t, TransactionStruct> gTransactionHistory;
static map<uint64_t,bool> gQueryResults;

static vector<uint32_t> gSchema;

//---------------------------------------------------------------------------
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
    //cerr << "Transaction: " << t.transactionId << endl;
    vector<TransOperation> operations;
    const char* reader=t.operations;

    // Delete all indicated tuples
    for (uint32_t index=0;index!=t.deleteCount;++index) {
        auto& o=*reinterpret_cast<const TransactionOperationDelete*>(reader);
        for (const uint64_t* key=o.keys,*keyLimit=key+o.rowCount;key!=keyLimit;++key) {
            auto& relation = gRelations[o.relationId];
            auto& rows = relation.insertedRows;
            auto lb = relation.insertedRows.find(*key);
            if (lb != rows.end()) {
                uint64_t rowid = lb->second[0];
                // update the relation transactions
                relation.transactions.push_back(TransStruct(t.transactionId, move(lb->second)));
                // update the transactions history record - so far we to not need tuple
                //operations.push_back(TransOperation(o.relationId, rowid, relation.transactions.back().tuple));
                operations.push_back(TransOperation(o.relationId, rowid));
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
            vector<uint64_t> tuple;
            tuple.insert(tuple.begin(),values,values+gSchema[o.relationId]);

            // add the tuple to this transaction operations and to the relations table
            // TODO - we might not have to store the tuple into the transaction history
            //operations.push_back(TransOperation(o.relationId, tuple[0], tuple));
            operations.push_back(TransOperation(o.relationId, tuple[0]));
            
            gRelations[o.relationId].transactions.push_back(TransStruct(t.transactionId, tuple));
            gRelations[o.relationId].insertedRows[values[0]]=move(tuple);
        }
        reader+=sizeof(TransactionOperationInsert)+(sizeof(uint64_t)*o.rowCount*gSchema[o.relationId]);
    }

    gTransactionHistory.insert(move(std::pair<uint64_t, TransactionStruct>(t.transactionId, TransactionStruct(move(operations))))); 
    //cerr << "Success Transaction: " << t.transactionId << endl;
}
//---------------------------------------------------------------------------
struct NullComb {
bool operator() (const TransStruct* p, std::nullptr_t target) {
    return p != target;
}
};
static void processValidationQueries(const ValidationQueries& v) {
    //cerr << "Validate: " << v.from << ":" << v.to << endl;
    TransStruct fromTRS(v.from);
    TransStruct toTRS(v.to);

    // try to put all the queries into a vector
    vector<LPQuery> queries;
    const char* qreader=v.queries;
    const Query *q;
    for (unsigned int i=0;i<v.queryCount;++i) {
        q=reinterpret_cast<const Query*>(qreader);
        queries.push_back(move(LPQuery(*q)));
        qreader+=sizeof(Query)+(sizeof(Query::Column)*q->columnCount);
    }

    // TODO - PROCESS the queries in order to take out common predicates on the same relations
    // TODO - to avoid multiple times the same query validation
    stable_sort(queries.begin(), queries.end());
    //for (unsigned int i=0;i<v.queryCount;++i) { cerr << queries[i].relationId << " "; } cerr << endl;
    queries.resize(std::distance(queries.begin(), unique(queries.begin(), queries.end())));
    //for (unsigned int i=0;i<v.queryCount;++i) { cerr << queries[i].relationId << " "; } cerr << endl;

    // TODO -  VERY NAIVE HERE - validate each query separately
    uint32_t lastRelId = gRelations.size()+1; // not valid relation id
    bool conflict=false;
    decltype(gRelations[0].transactions)::iterator transFrom;
    decltype(gRelations[0].transactions)::iterator transTo;
    for (unsigned int index=0,qsz=queries.size();index<qsz;++index) {
        auto& q=queries[index];
        // avoid searching for the range of transactions too many times 
        if (q.relationId != lastRelId) {
            lastRelId = q.relationId;
            auto& relation = gRelations[q.relationId];
            auto& transactions = relation.transactions;
            transFrom = std::lower_bound(transactions.begin(), transactions.end(), fromTRS, TRSLessThan);
            transTo = std::upper_bound(transactions.begin(), transactions.end(), toTRS, TRSLessThan);
        }
        
        for(auto iter=transFrom;iter!=transTo;++iter) {
            auto& tuple = iter->tuple;
            bool match=true;
            for (auto c=q.columns.begin(),cLimit=q.columns.end();c!=cLimit;++c) {
                // make the actual check
                uint64_t tupleValue = tuple[c->column];
                uint64_t queryValue=c->value;
                bool result=false;
                    switch (c->op) {
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
            }
            if (match) {
                conflict=true;
                break;
            }    
        } // end of all transactions for this relation query
        if (conflict) break;
    }
    gQueryResults[v.validationId]=conflict;
    //cerr << "Success Validate: " << v.validationId << " ::" << v.from << ":" << v.to << " result: " << conflict << endl;
}
//---------------------------------------------------------------------------
static void processFlush(const Flush& f) {
    // TODO - NEEDS A COMPLETE REFACTORING
    while ((!gQueryResults.empty())&&((*gQueryResults.begin()).first<=f.validationId)) {
        char c='0'+(*gQueryResults.begin()).second;
        cout.write(&c,1);
        gQueryResults.erase(gQueryResults.begin());
    }
    cout.flush();
}
//---------------------------------------------------------------------------

static void processForget(const Forget& f) {
    //cerr << "Forget: " << f.transactionId << endl;
    //auto start = std::chrono::system_clock::now();
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
    //auto end = std::chrono::system_clock::now();
    //cerr << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << endl;
}


//---------------------------------------------------------------------------
// Read the message body and cast it to the desired type
template<typename Type> static const Type& readBody(istream& in,vector<char>& buffer,uint32_t len) {
    buffer.resize(len);
    in.read(buffer.data(),len);
    return *reinterpret_cast<const Type*>(buffer.data());
}
//---------------------------------------------------------------------------
int main()
{
    vector<char> message;
    while (true) {
        // Retrieve the message
        MessageHead head;
        cin.read(reinterpret_cast<char*>(&head),sizeof(head));
        if (!cin) { cerr << "read error" << endl; abort(); } // crude error handling, should never happen

        // And interpret it
        switch (head.type) {
            case MessageHead::Transaction: processTransaction(readBody<Transaction>(cin,message,head.messageLen)); break;
            case MessageHead::ValidationQueries: processValidationQueries(readBody<ValidationQueries>(cin,message,head.messageLen)); break;
            case MessageHead::Flush: processFlush(readBody<Flush>(cin,message,head.messageLen)); break;
            case MessageHead::Forget: processForget(readBody<Forget>(cin,message,head.messageLen)); break;
            case MessageHead::Done: return 0;
            case MessageHead::DefineSchema: processDefineSchema(readBody<DefineSchema>(cin,message,head.messageLen)); break;
            default: cerr << "malformed message" << endl; abort(); // crude error handling, should never happen
        }
    }
}
//---------------------------------------------------------------------------
