//////////////////////////////////////////////////////
ACM SIGMOD 2015 Programming Contest
http://db.in.tum.de/sigmod15contest/index.html

(1) Team name: CStrings

(2) Members:

* Lambros Petrou | lambrospetrou@gmail.com, University of Oxford, UK, MSc Computer Science
* George Koumettou | gkoume01@gmail.com, Royal Holloway University of London (RHUL), UK, MSc Information Security
* Marios Mintzis | mariosmintzis@hotmail.com, University College London, UK, MSc Web Science and Big Data Analytics

All 3 members graduated in 2014 from University of Cyprus (UCY) with a BSc of Computer Science.

(3) Supervisor name: (-)

(4) Summary of techniques
*) Dedicated thread did the reading so the main thread always had a queued message to parse
*) Main thread did the message parsing:
    validation messages: just queued without processing
    transaction messages: parsed only to distribute part of the transaction to the relations affected
    flush/forget messages: block -> initiate execution of all the transactions and then all the validations queued ->
                            case forget -> delete part of the index
                            case flush -> print out results

*) Main execution consisted of 4 stages
    1) process transactions on each relation concurrently
    2) for each relation finished step 1, we process each column concurrently to build our Column-based Index
    3) concurrently parse all validation messages and validate the queries -> prune out the invalid queries
        - sort query predicates by operator and then distribute them to the relation/column 
            of their 1st predicate (after sorting) to build the query index
        - queries without an Equality (==) operator were left in the validations queue and not in the query index
    4) all threads concurrently take columns from the query index and execute all the queries of that column 
        using the column index created at step 2.
    5) process any remaining queries from the original validation queue (queries without any == operator)

*) The key features of our solution were:
    1) the fact that we executed all the queries related to each column together so we were using the same column index for thousands of queries. And since the majority of the queries had at least 1 == operator the ones left out for further processing (step 5) were minimal and were executed only if the validation has not been found conflicting yet.
    
    2) We designed a custom Column-based Index for the transactions and the tuples they deleted/inserted. Each column had buckets of transactions and the buckets were sorted by the transaction Ids they contained in order to allow search only in the transaction range specified for each validation. Additionally, each bucket was filling up until either:
        a) the number of transactions contained passed our threshold
        b) the number of tuples contained passed our threshold
    and as a result we always had balanced buckets with just enough tuples to being fast to process. Each bucket inside had the tuples sorted by their value in that specific column, therefore with just a single binary search we could get all the tuples that had the value we want.
    
Open source libraries:
- Subroutine library by Agner: http://www.agner.org/optimize/
- C++ B-tree by Google: https://code.google.com/p/cpp-btree/
- Concurrent Queue: https://github.com/cameron314/readerwriterqueue

Most important thing - the contest was really fun and we learnt a ton from our mistakes!!!

//////////////////////////////////////////////////////
Compile the program with
./compile.sh

Test the program with ( this will only test the small basic.test file)
./test.sh

For uploading your submission you can use the
./package.sh
which will bundle the following:
   src/
   run.sh
   test.sh

Before packaging please make sure that you remove all binaries
from the src/ folder
