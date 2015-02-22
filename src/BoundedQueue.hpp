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

    BoundedQueue(uint64_t sz) : mMaxSize(sz), InvalidPtr(sz+1), mlUnused(InvalidPtr), mlAvailable(InvalidPtr), 
        mlReqEnq(InvalidPtr), mlReqDeq(InvalidPtr), mNodes(std::vector<Node>(mMaxSize)) {
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
    }

    BQResult reqNextEnq() { 
        std::unique_lock<std::mutex> lk(mMutex);
        debugInfo("req-enq");
        mCondFull.wait(lk, [this]{return mlUnused != InvalidPtr;});
        // transfer the node from Unused list to ReqEnq list
        NodePtr cN = disconnectHead(mlUnused);
        // connect node to mReqEnq
        setHead(cN, mlReqEnq);
        lk.unlock();
        return BQResult(cN, &mNodes[cN].value); 
    }

    void registerEnq(uint64_t id) {
        // might not be necessary
        std::lock_guard<std::mutex> lk(mMutex);
        debugInfo("register-enq");
        // transfer node from mReqEnq to Available
        NodePtr nN = mNodes[id].refId;
        // remove it from ReqEnq
        disconnectNode(nN, mlReqEnq);
        // connect it to Available
        setHead(nN, mlAvailable);
        mCondEmpty.notify_one();
    }
    
    BQResult reqNextDeq() {
        std::unique_lock<std::mutex> lk(mMutex);
        debugInfo("req-deq");
        mCondEmpty.wait(lk, [this]{return mlAvailable != InvalidPtr;});
        // transfer the node from Available list to ReqDeq list
        NodePtr cN = disconnectHead(mlAvailable);
        // connect node to mReqDeq
        setHead(cN, mlReqDeq);
        lk.unlock();
        return BQResult(cN, &mNodes[cN].value);
    }
    void registerDeq(uint64_t id) {
        std::lock_guard<std::mutex> lk(mMutex);
        debugInfo("register-deq");
        // transfer node from mReqDeq to Unused
        NodePtr cN = mNodes[id].refId;
        // remove it from ReqEnq
        disconnectNode(cN, mlReqDeq);
        // connect it to Available
        setHead(cN, mlUnused);
        mCondFull.notify_one();    
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

    void disconnectNode(NodePtr cN, NodePtr& l) {
        if (cN == l) {
            // we are removing the head
            l = mNodes[cN].next;
            if (l != InvalidPtr) mNodes[l].prev = InvalidPtr;
        } else { 
            if (mNodes[cN].next != InvalidPtr) mNodes[mNodes[cN].next].prev = mNodes[cN].prev;
            if (mNodes[cN].prev != InvalidPtr) mNodes[mNodes[cN].prev].next = mNodes[cN].next;
        }
    }

    NodePtr disconnectHead(NodePtr &l) {
        if (l == InvalidPtr) return InvalidPtr;
        NodePtr t = l;
        l = mNodes[l].next;
        return t;
    }

    void setHead(NodePtr nN, NodePtr& l) { 
        mNodes[nN].next = l;
        mNodes[nN].prev = InvalidPtr;
        if (l != InvalidPtr) mNodes[l].prev = nN;
        l = nN;
    }  
};

#endif
