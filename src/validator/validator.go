// SIGMOD Programming Contest 2015 - stub implementation
//
// This code is intended to illustrate the given challenge, it
// is not particular efficient. It is not guaranteed to be correct!
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

package main

import (
   "bufio"
   "encoding/binary"
   "fmt"
   "io"
   "os"
)

type Row struct {
   transaction uint64
   columns []uint64
}

type Relation struct {
   columns uint32
   inserted map[uint64]Row
   deleted []Row
}

type Validation struct {
   id uint64
   success bool
}

type Predicate struct {
   column uint32
   operation uint32
   reference uint64
}

type Validator struct {
   relations []Relation
   validations []Validation
   last_forget uint64
}

func (v *Validator) DefineSchema(r io.Reader) {
   var count uint32
   _ = binary.Read(r, binary.LittleEndian, &count)
   for i := 0; i < int(count); i++ {
      var relation Relation
      relation.inserted = make(map[uint64]Row)
      _ = binary.Read(r, binary.LittleEndian, &relation.columns)
      v.relations = append(v.relations, relation)
   }
}

func (v *Validator) Transaction(r io.Reader) {
   var transaction uint64
   var deletions uint32
   var insertions uint32
   _ = binary.Read(r, binary.LittleEndian, &transaction)
   _ = binary.Read(r, binary.LittleEndian, &deletions)
   _ = binary.Read(r, binary.LittleEndian, &insertions)

   for i := 0; i < int(deletions); i++ {
      var relation uint32
      var rows uint32
      _ = binary.Read(r, binary.LittleEndian, &relation)
      _ = binary.Read(r, binary.LittleEndian, &rows)

      for j := 0; j < int(rows); j++ {
         var key uint64
         _ = binary.Read(r, binary.LittleEndian, &key)
         v.DeleteRow(transaction, relation, key)
      }
   }

   for i := 0; i < int(insertions); i++ {
      var relation uint32
      var rows uint32
      _ = binary.Read(r, binary.LittleEndian, &relation)
      _ = binary.Read(r, binary.LittleEndian, &rows)

      for j := 0; j < int(rows); j++ {
         v.InsertRow(transaction, relation, r)
      }
   }
}

func (v *Validator) Validate(r io.Reader) {
   var validation Validation
   var validation_id uint64
   var from_tx uint64
   var to_tx uint64
   var queries uint32
   _ = binary.Read(r, binary.LittleEndian, &validation_id)
   _ = binary.Read(r, binary.LittleEndian, &from_tx)
   _ = binary.Read(r, binary.LittleEndian, &to_tx)
   _ = binary.Read(r, binary.LittleEndian, &queries)


   validation.id = validation_id
   validation.success = true

   for i := 0; i < int(queries); i++ {
      var predicates []Predicate
      var relation uint32
      var count uint32
      _ = binary.Read(r, binary.LittleEndian, &relation)
      _ = binary.Read(r, binary.LittleEndian, &count)
      if count == 0 {
         var predicate Predicate
         predicate.column = 0
         predicate.operation = 6
         predicate.reference = 0
         predicates = append(predicates, predicate)
      } else {
         for j := 0; j < int(count); j++ {
            var predicate Predicate
            _ = binary.Read(r, binary.LittleEndian, &predicate.column)
            _ = binary.Read(r, binary.LittleEndian, &predicate.operation)
            _ = binary.Read(r, binary.LittleEndian, &predicate.reference)
            predicates = append(predicates, predicate)
         }
      }

      if validation.success {
         result := v.ValidatePredicates(from_tx, to_tx, relation, predicates)
         if !result {
            validation.success = false
         }
      }
   }


   v.validations = append(v.validations, validation)
}

func (v *Validator) ValidatePredicates(from_tx uint64, to_tx uint64, relation uint32, predicates []Predicate) bool {
   for _, row := range v.relations[relation].inserted {
      if row.transaction >= from_tx && row.transaction <= to_tx {
         all_match := true
         for _, predicate := range predicates {
            result := ValidatePredicate(row, predicate)
            if !result {
               all_match = false
               break
            }
         }
         if all_match {
            return false
         }
      }
   }

   for _, row := range v.relations[relation].deleted {
      if row.transaction >= from_tx && row.transaction <= to_tx {
         all_match := true
         for _, predicate := range predicates {
            result := ValidatePredicate(row, predicate)
            if !result {
               all_match = false
               break
            }
         }
         if all_match {
            return false
         }
      }
   }

   return true
}

func ValidatePredicate(row Row, predicate Predicate) bool {
   column := int(predicate.column)
   switch predicate.operation {
      case 0: return row.columns[column]==predicate.reference
      case 1: return row.columns[column]!=predicate.reference
      case 2: return row.columns[column]<predicate.reference
      case 3: return row.columns[column]<=predicate.reference
      case 4: return row.columns[column]>predicate.reference
      case 5: return row.columns[column]>=predicate.reference
      case 6: return true
      default: panic("Invalid operator in predicate")
   }

   return false
}

func (v *Validator) Flush(r io.Reader) {
      var flush_id uint64
      _ = binary.Read(r, binary.LittleEndian, &flush_id)


      fmt.Fprintf(os.Stderr, "Flush up to %d\n", flush_id)

      var remaining []Validation
      for _, validation := range v.validations {
         if validation.id <= flush_id {
            if validation.success {
               os.Stdout.Write([]byte{byte('0')})
            } else {
               os.Stdout.Write([]byte{byte('1')})
            }
         } else {
            remaining = append(remaining, validation)
         }
      }
      v.validations = remaining
}

func (v *Validator) Forget(r io.Reader) {
   var forget_id uint64
   _ = binary.Read(r, binary.LittleEndian, &forget_id)
   
   fmt.Fprintf(os.Stderr, "\nForget up to: %d\n", forget_id)
   // TODO: Implement
}

func (v *Validator) DeleteRow(transaction uint64, relation uint32, key uint64) {
   row, ok := v.relations[relation].inserted[key]
   if ok {
      v.relations[relation].deleted = append(v.relations[relation].deleted, row)
      row.transaction = transaction
      v.relations[relation].deleted = append(v.relations[relation].deleted, row)
      delete(v.relations[relation].inserted, key)
   }
}

func (v *Validator) InsertRow(transaction uint64, relation uint32, r io.Reader) {
   var row Row
   row.transaction = transaction
   for i := 0; i < int(v.relations[relation].columns); i++ {
      var value uint64
      _ = binary.Read(r, binary.LittleEndian, &value)
      row.columns = append(row.columns, value)
   }
   v.relations[relation].inserted[row.columns[0]] = row
}

func read_requests(r io.Reader) {
   var validator Validator

   transactions := 0
   validations := 0

   for {
      var rlen uint32
      var rtype uint32

      err := binary.Read(r, binary.LittleEndian, &rlen)
      if err != nil {
         panic("Error during reading")
      }
      _ = binary.Read(r, binary.LittleEndian, &rtype)

      switch rtype {
         case 0 : {fmt.Fprintf(os.Stderr, "\nTransactions: %d\nValidations: %d\n", transactions, validations); return }
         case 1 : validator.DefineSchema(r)
         case 2 : { validator.Transaction(r); transactions += 1 }
         case 3 : { validator.Validate(r); validations += 1 }
         case 4 : validator.Flush(r)
         case 5 : { fmt.Fprintf(os.Stderr, "\nTransactions so far: %d", transactions); validator.Forget(r) }
         default : panic("Invalid message type")
      }
   }
   
}
func main() {
    reader := bufio.NewReader(os.Stdin)
    read_requests(reader)
}
