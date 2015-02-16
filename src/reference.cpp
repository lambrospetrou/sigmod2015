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
#include <vector>
#include <cassert>
#include <cstdint>
#include <algorithm>
#include <utility>
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

// Custom data structures to hold data
struct CTransStruct{
    uint64_t trans_id;
    uint64_t value;
    CTransStruct(uint64_t tid, uint64_t v): trans_id(tid), value(v) {}
};
struct CTRSLessThan_t {
    bool operator() (const CTransStruct& left, const CTransStruct& right) {
        return left.trans_id < right.trans_id;
    }
    bool operator() (const CTransStruct& o, uint64_t target) {
        return o.trans_id < target;
    }
    bool operator() (uint64_t target, const CTransStruct& o) {
        return target < o.trans_id;
    }
} CTRSLessThan;

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
    map<uint32_t, vector<uint64_t>> insertedRows;
};
static vector<RelationStruct> gRelations;

struct TransOperation {
    uint32_t rel_id;
    uint64_t row_id;
    vector<uint64_t> tuple;
    TransOperation(uint64_t relid, uint32_t rowid, vector<uint64_t> t):
        rel_id(relid), row_id(rowid), tuple(t) {}
};
struct TransactionStruct {
    //    uint64_t trans_id; // we will use direct indexing
    vector<TransOperation> operations;
    TransactionStruct(vector<TransOperation> ops) : operations(ops) {}
};
static map<uint64_t, TransactionStruct> gTransactionHistory;
static map<uint64_t,bool> gQueryResults;

static vector<uint32_t> gSchema;

/*
   static vector<uint32_t> schema;
   static vector<map<uint32_t,vector<uint64_t>>> relations;
//---------------------------------------------------------------------------
static map<uint64_t,vector<pair<uint32_t,vector<uint64_t>>>> transactionHistory;
static map<uint64_t,bool> queryResults;
 */
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

    //vector<pair<uint32_t,vector<uint64_t>>> operations;
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
                // update the transactions history record
                operations.push_back(TransOperation(o.relationId, rowid, relation.transactions.back().tuple));
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
            operations.push_back(TransOperation(o.relationId, tuple[0], tuple));
            
            gRelations[o.relationId].transactions.push_back(TransStruct(t.transactionId, tuple));
            gRelations[o.relationId].insertedRows[values[0]]=move(tuple);
        }
        reader+=sizeof(TransactionOperationInsert)+(sizeof(uint64_t)*o.rowCount*gSchema[o.relationId]);
    }

    gTransactionHistory.insert(std::pair<uint64_t, TransactionStruct>(t.transactionId, TransactionStruct(move(operations)))); 
    //cerr << "Success Transaction: " << t.transactionId << endl;
}
//---------------------------------------------------------------------------
static void processValidationQueries(const ValidationQueries& v) {
    //cerr << "Validate: " << v.from << ":" << v.to << endl;
    TransStruct fromTRS(v.from);
    TransStruct toTRS(v.to);

    // TODO - VERY NAIVE HERE - validate each query separately
    bool conflict=false;
    const char* reader=v.queries;
    for (unsigned int index=0;index<v.queryCount;++index) {
        auto& q=*reinterpret_cast<const Query*>(reader);

        auto& relation = gRelations[q.relationId];
        auto& transactions = relation.transactions;
        auto transFrom = std::lower_bound(transactions.begin(), transactions.end(), fromTRS, TRSLessThan);
        auto transTo = std::upper_bound(transactions.begin(), transactions.end(), toTRS, TRSLessThan);

        for(;transFrom!=transTo;++transFrom) {
            auto& tuple = transFrom->tuple;
            bool match=true;
            for (auto c=q.columns,cLimit=c+q.columnCount;c!=cLimit;++c) {
                //cerr << "\nquery[" << index << "] rel[" << q.relationId << "] col[" << c->column << "] op[" << c->op << "] qv[" << c->value << "]" << endl;

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
        reader+=sizeof(Query)+(sizeof(Query::Column)*q.columnCount);
    }
    gQueryResults[v.validationId]=conflict;
    //cerr << "Success Validate: " << v.validationId << " ::" << v.from << ":" << v.to << " result: " << conflict << endl;
    //if (v.validationId == 21) exit(1);
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
    auto iend = gTransactionHistory.upper_bound(f.transactionId);
    TransStruct fstruct(f.transactionId);

    // Delete the transactions from inside the columns in the relations
    vector<bool> relationsCleaned(gRelations.size(), false);
    for(auto ctr=gTransactionHistory.begin(); ctr!=iend; ++ctr) {
        //cerr << ":: each transaction: " << ctr->second.operations.size() << endl;
        for(auto cop=ctr->second.operations.begin(), opend=ctr->second.operations.end(); cop!=opend; ++cop ) {
            //cerr << ":: each operation: " << cop->rel_id << ":" << gRelations.size() << endl;
            if (!relationsCleaned[cop->rel_id]) {
                relationsCleaned[cop->rel_id] = true;
                // delete this transaction from the lastRel columns
                auto& transactions = gRelations[cop->rel_id].transactions;
                transactions.erase(transactions.begin(), 
                    upper_bound(transactions.begin(), transactions.end(), fstruct, TRSLessThan));
            }
        }
    }

    // then delete the transactions from the transaction history
    gTransactionHistory.erase(gTransactionHistory.begin(), iend);
    //cerr << ":: Success Forget: " << f.transactionId << endl;
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
            case MessageHead::Done: return 0;
            case MessageHead::DefineSchema: processDefineSchema(readBody<DefineSchema>(cin,message,head.messageLen)); break;
            case MessageHead::Transaction: processTransaction(readBody<Transaction>(cin,message,head.messageLen)); break;
            case MessageHead::ValidationQueries: processValidationQueries(readBody<ValidationQueries>(cin,message,head.messageLen)); break;
            case MessageHead::Flush: processFlush(readBody<Flush>(cin,message,head.messageLen)); break;
            case MessageHead::Forget: processForget(readBody<Forget>(cin,message,head.messageLen)); break;
            default: cerr << "malformed message" << endl; abort(); // crude error handling, should never happen
        }
    }
}
//---------------------------------------------------------------------------
