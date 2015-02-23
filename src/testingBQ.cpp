#include <thread>
#include <iostream>
#include <functional>
#include <cstdint>
#include <chrono>

#include "BoundedQueue.hpp"
#include "MultiTaskPool.hpp"

using namespace std;


void producer2(BoundedQueue<uint64_t>* Q) {
    uint64_t num = 0;
    while(true) {
        auto res = Q->reqNextEnq();
        cout << "enq id: " << res.refId << endl;
        *res.value = num++;
        Q->registerEnq(res.refId);
    }
}

void producer(BoundedQueue<char>* Q) {
    while(true) {
        auto res = Q->reqNextEnq();
        cout << "enq id: " << res.refId << endl;
        *res.value = 'a';
        Q->registerEnq(res.refId);
        
        auto res1 = Q->reqNextEnq();
        cout << "enq id: " << res.refId << endl;
        auto res2 = Q->reqNextEnq();
        cout << "enq id: " << res.refId << endl;
        auto res3 = Q->reqNextEnq();
        cout << "enq id: " << res.refId << endl;

        *res1.value = 'b';
        *res2.value = 'c';
        *res3.value = 'd';

        Q->registerEnq(res1.refId);
        Q->registerEnq(res2.refId);
        Q->registerEnq(res3.refId);
    }
}

void consumer(BoundedQueue<uint64_t>* Q) {
    while(true) {
        auto res = Q->reqNextDeq();
        //auto res1 = Q->reqNextDeq();
        //auto res2 = Q->reqNextDeq();
        cout << "deq id: " << res.refId << " deq val: " << *res.value << endl;
        this_thread::sleep_for(chrono::milliseconds(20));
        Q->registerDeq(res.refId);
        //Q->registerDeq(res1.refId);
        //Q->registerDeq(res2.refId);
    }
}   

void workerFunc(uint32_t nThreads, uint32_t tid, void *args) {
    int* iptr = static_cast<int*>(args);
    cerr << "worker " << tid << "/" << nThreads << " args: " << *iptr << endl;
    delete iptr;
}

int main() {
/*
    BoundedQueue<uint64_t> Q(100);

    thread t1(producer2, &Q);
    thread t2(consumer, &Q);
    t1.join();
    t2.join();
*/

    MultiTaskPool pool(4);
    pool.initThreads();
    pool.startAll();
    for(int i=0; i<100; i++) {
        cerr << "job added" << pool.addTask(workerFunc, static_cast<void*>(new int(i))) << endl;
    }
    pool.waitAllAndStop();
    cerr << "SHOULD SYNC" << endl;
    for(int i=0; i<100; i++) {
        cerr << "job added" << pool.addTask(workerFunc, static_cast<void*>(new int(i))) << endl;
    }
    pool.startAll();
    pool.waitAll();
    pool.destroy();

    return 0;
}









