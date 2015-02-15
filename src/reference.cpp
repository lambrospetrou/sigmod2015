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
static vector<uint32_t> schema;
static vector<map<uint32_t,vector<uint64_t>>> relations;
//---------------------------------------------------------------------------
static map<uint64_t,vector<pair<uint32_t,vector<uint64_t>>>> transactionHistory;
static map<uint64_t,bool> queryResults;
//---------------------------------------------------------------------------
static void processDefineSchema(const DefineSchema& d)
{
   schema.clear();
   schema.insert(schema.begin(),d.columnCounts,d.columnCounts+d.relationCount);
   relations.clear();
   relations.resize(d.relationCount);
}
//---------------------------------------------------------------------------
static void processTransaction(const Transaction& t)
{
   vector<pair<uint32_t,vector<uint64_t>>> operations;
   const char* reader=t.operations;

   // Delete all indicated tuples
   for (uint32_t index=0;index!=t.deleteCount;++index) {
      auto& o=*reinterpret_cast<const TransactionOperationDelete*>(reader);
      for (const uint64_t* key=o.keys,*keyLimit=key+o.rowCount;key!=keyLimit;++key) {
         if (relations[o.relationId].count(*key)) {
            operations.push_back(pair<uint32_t,vector<uint64_t>>(o.relationId,move(relations[o.relationId][*key])));
            relations[o.relationId].erase(*key);
         }
      }
      reader+=sizeof(TransactionOperationDelete)+(sizeof(uint64_t)*o.rowCount);
   }

   // Insert new tuples
   for (uint32_t index=0;index!=t.insertCount;++index) {
      auto& o=*reinterpret_cast<const TransactionOperationInsert*>(reader);
      for (const uint64_t* values=o.values,*valuesLimit=values+(o.rowCount*schema[o.relationId]);values!=valuesLimit;values+=schema[o.relationId]) {
         vector<uint64_t> tuple;
         tuple.insert(tuple.begin(),values,values+schema[o.relationId]);
         operations.push_back(pair<uint32_t,vector<uint64_t>>(o.relationId,tuple));
         relations[o.relationId][values[0]]=move(tuple);
      }
      reader+=sizeof(TransactionOperationInsert)+(sizeof(uint64_t)*o.rowCount*schema[o.relationId]);
   }

   transactionHistory[t.transactionId]=move(operations);
}
//---------------------------------------------------------------------------
static void processValidationQueries(const ValidationQueries& v)
{
   auto from=transactionHistory.lower_bound(v.from);
   auto to=transactionHistory.upper_bound(v.to);
   bool conflict=false;
   const char* reader=v.queries;
   for (unsigned index=0;index!=v.queryCount;++index) {
      auto& q=*reinterpret_cast<const Query*>(reader);
      for (auto iter=from;iter!=to;++iter) {
         for (auto& op:(*iter).second) {
            // Check if the relation is the same
            if (op.first!=q.relationId)
               continue;

            // Check if all predicates are satisfied
            auto& tuple=op.second;
            bool match=true;
            for (auto c=q.columns,cLimit=c+q.columnCount;c!=cLimit;++c) {
               uint64_t tupleValue=tuple[c->column],queryValue=c->value;
               bool result=false;
               switch (c->op) {
                  case Query::Column::Equal: result=(tupleValue==queryValue); break;
                  case Query::Column::NotEqual: result=(tupleValue!=queryValue); break;
                  case Query::Column::Less: result=(tupleValue<queryValue); break;
                  case Query::Column::LessOrEqual: result=(tupleValue<=queryValue); break;
                  case Query::Column::Greater: result=(tupleValue>queryValue); break;
                  case Query::Column::GreaterOrEqual: result=(tupleValue>=queryValue); break;
               }
               if (!result) { match=false; break; }
            }
            if (match) {
               conflict=true;
               break;
            }
         }
      }
      reader+=sizeof(Query)+(sizeof(Query::Column)*q.columnCount);
   }

   queryResults[v.validationId]=conflict;
}
//---------------------------------------------------------------------------
static void processFlush(const Flush& f)
{
   while ((!queryResults.empty())&&((*queryResults.begin()).first<=f.validationId)) {
      char c='0'+(*queryResults.begin()).second;
      cout.write(&c,1);
      queryResults.erase(queryResults.begin());
   }
   cout.flush();
}
//---------------------------------------------------------------------------
static void processForget(const Forget& f)
{
   while ((!transactionHistory.empty())&&((*transactionHistory.begin()).first<=f.transactionId))
      transactionHistory.erase(transactionHistory.begin());
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
