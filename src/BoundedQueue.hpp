#ifndef __LP_BOUNDEDQUEUE__
#define __LP_BOUNDEDQUEUE__

#include <vector>
#include <mutex>
#include <condition_variable>
#include <cstdint>
#include <chrono>
#include <thread>
#include <functional>
#include <memory>
#include <iostream>
#include <cstdint>
#include <string>

template<class T>
class BoundedQueue {
private: 
    typedef uint64_t NodePtr;

    struct Node {
        NodePtr next;
        NodePtr prev;
        NodePtr refId;
        T value;
    };

public:

    struct BQResult {
        uint64_t refId;
        T* value;
        BQResult(uint64_t id, T* v) : refId(id), value(v){}
    };

    BoundedQueue(uint64_t sz) : mMaxSize(sz), InvalidPtr(sz+1), 
    mlUnused(InvalidPtr), mlAvailable(InvalidPtr), mlReqEnq(InvalidPtr), mlReqDeq(InvalidPtr), 
    mlUnusedTail(InvalidPtr), mlAvailableTail(InvalidPtr), mlReqEnqTail(InvalidPtr), mlReqDeqTail(InvalidPtr), 
    mNodes(std::vector<Node>(mMaxSize)) {
            if (sz == 0) return;
            mNodes[0].prev = InvalidPtr;
            mNodes[0].refId = 0;
            for (uint64_t i=1; i<sz; ++i) {
                mNodes[i].prev = i-1;
                mNodes[i-1].next = i;
                mNodes[i].refId = i;
            }
            mNodes[sz-1].next = InvalidPtr;
            mlUnused = 0;
            mlUnusedTail = sz-1;
            mReservedDeq = 0;
            mCurEnqueued = 0;
    }

    BQResult reqNextEnq() { 
        std::unique_lock<std::mutex> lk(mMutex);
        //debugInfo("req-enq");
        mCondFull.wait(lk, [this]{return mlUnused != InvalidPtr;});
        // transfer the node from Unused list to ReqEnq list
        NodePtr cN = disconnectHead(mlUnused, mlUnusedTail);
        // connect node to mReqEnq
        append(cN, mlReqEnq, mlReqEnqTail);
        lk.unlock();
        return BQResult(cN, &mNodes[cN].value); 
    }

    void registerEnq(uint64_t id) {
        // might not be necessary
        std::lock_guard<std::mutex> lk(mMutex);
        //debugInfo("register-enq");
        // transfer node from mReqEnq to Available
        NodePtr nN = mNodes[id].refId;
        // remove it from ReqEnq
        disconnectNode(nN, mlReqEnq, mlReqEnqTail);
        // connect it to Available
        append(nN, mlAvailable, mlAvailableTail);
        ++mCurEnqueued;
        mCondEmpty.notify_one();
    }
    
    BQResult reqNextDeq() {
        std::unique_lock<std::mutex> lk(mMutex);
        ++mReservedDeq;
        //debugInfo("req-deq");
        mCondEmpty.wait(lk, [this]{return mlAvailable != InvalidPtr;});
        // transfer the node from Available list to ReqDeq list
        NodePtr cN = disconnectHead(mlAvailable, mlAvailableTail);
        // connect node to mReqDeq
        append(cN, mlReqDeq, mlReqDeqTail);
        lk.unlock();
        return BQResult(cN, &mNodes[cN].value);
    }
    void registerDeq(uint64_t id) {
        std::lock_guard<std::mutex> lk(mMutex);
        --mReservedDeq;
        //debugInfo("register-deq");
        // transfer node from mReqDeq to Unused
        NodePtr cN = mNodes[id].refId;
        // remove it from ReqEnq
        disconnectNode(cN, mlReqDeq, mlReqDeqTail);
        // connect it to Available
        append(cN, mlUnused, mlUnusedTail);
        if (--mCurEnqueued <= (mMaxSize>>1)) mCondFull.notify_one();    
    }
  
    void debugInfo(const std::string& s) {
        std::cerr << std::endl << s << std::endl;
        std::cerr << "Unused:"; cerrList(mlUnused);
        std::cerr << "Avail:";  cerrList(mlAvailable);
        std::cerr << "ReqEnq:"; cerrList(mlReqEnq);
        std::cerr << "ReqDeq:"; cerrList(mlReqDeq);
    }

private:

    const uint64_t mMaxSize;
    const uint64_t InvalidPtr; 
    
    NodePtr mlUnused;
    NodePtr mlAvailable;
    NodePtr mlReqEnq;
    NodePtr mlReqDeq;
    
    NodePtr mlUnusedTail;
    NodePtr mlAvailableTail;
    NodePtr mlReqEnqTail;
    NodePtr mlReqDeqTail;
   
    uint64_t mReservedDeq;
    uint64_t mCurEnqueued;

    std::vector<Node> mNodes;
    std::mutex mMutex;
    std::condition_variable mCondFull;
    std::condition_variable mCondEmpty;

    void cerrList(NodePtr n) {
        while (n != InvalidPtr) {
            std::cerr << "[" << n << "] ";
            n = mNodes[n].next;
        }
        std::cerr << std::endl;
    }

    void disconnectNode(NodePtr cN, NodePtr& lH, NodePtr& lT) {
        if (cN == InvalidPtr) return;
        if (cN == lH) {
            disconnectHead(lH, lT);
        } else if (cN == lT) {
            disconnectTail(lH, lT);  
        } else { 
            if (mNodes[cN].next != InvalidPtr) mNodes[mNodes[cN].next].prev = mNodes[cN].prev;
            if (mNodes[cN].prev != InvalidPtr) mNodes[mNodes[cN].prev].next = mNodes[cN].next;
        }
    }

    NodePtr disconnectHead(NodePtr &lH, NodePtr& lT) {
        if (lH == InvalidPtr) return InvalidPtr;
        NodePtr t = lH;
        lH = mNodes[lH].next;
        if (lH == InvalidPtr) lT = InvalidPtr;
        else mNodes[lH].prev = InvalidPtr;
        return t;
    }
    NodePtr disconnectTail(NodePtr &lH, NodePtr& lT) {
        if (lT == InvalidPtr) return InvalidPtr;
        NodePtr t = lT;
        lT = mNodes[lT].prev;
        if (lT == InvalidPtr) lH = InvalidPtr;
        else mNodes[lT].next = InvalidPtr;
        return t;
    }

    void prepend(NodePtr nN, NodePtr& lH, NodePtr& lT) { 
        if (lH == InvalidPtr) {
            // empty list
            mNodes[nN].prev = InvalidPtr;
            mNodes[nN].next = InvalidPtr;
            lH = nN; lT = nN;
        } else {
            // there is a head already
            mNodes[nN].next = lH;
            mNodes[nN].prev = InvalidPtr;
            mNodes[lH].prev = nN;
            lH = nN;
        }
    }  
    void append(NodePtr nN, NodePtr& lH, NodePtr& lT) { 
        if (lH == InvalidPtr) {
            // empty list so set both to this node
            lH = nN; lT = nN;
            mNodes[nN].prev = InvalidPtr;
            mNodes[nN].next = InvalidPtr;
        } else {
            // the list is not empty
            mNodes[nN].prev = lT;
            mNodes[nN].next = InvalidPtr;
            mNodes[lT].next = nN;
            lT = nN;
        }
    }
};

#endif
