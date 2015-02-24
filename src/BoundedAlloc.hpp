#ifndef __LP_BOUNDEDALLOC__
#define __LP_BOUNDEDALLOC__

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
#include <utility>

template<class T>
class BoundedAlloc {
    public:
        struct BAResult {
            uint64_t refId;
            T* value;
            BAResult(uint64_t id, T* v) : refId(id), value(v){}
            BAResult() : refId(UINT64_MAX), value(nullptr){}
        };
    private: 

        struct Node {
            Node* next;
            Node* prev;
            uint64_t refId;
            T value;

            BAResult result;
        };
        //typedef uint64_t NodePtr;
        typedef Node* NodePtr;

    public:

        BoundedAlloc(uint64_t sz) : mMaxSize(sz), InvalidPtr(nullptr), 
        mlUnused(InvalidPtr), mlReserved(InvalidPtr), mlUnusedTail(InvalidPtr), mlReservedTail(InvalidPtr), 
        mNodes(std::vector<Node>(mMaxSize)) {
            if (sz == 0) return;
            mNodes[0].result.refId = 0;
            mNodes[0].result.value = &mNodes[0].value;
            mNodes[0].prev = InvalidPtr;
            mNodes[0].refId = 0;
            for (uint64_t i=1; i<sz; ++i) {
                mNodes[i].prev = &mNodes[i-1];
                mNodes[i-1].next = &mNodes[i];
                mNodes[i].refId = i;
                mNodes[i].result.refId = i;
                mNodes[i].result.value = &mNodes[i].value;
            }
            mNodes[sz-1].next = InvalidPtr;
            mNodes[sz-1].result.refId = sz-1;
            mNodes[sz-1].result.value = &mNodes[sz-1].value;
            mlUnused = &mNodes[0];
            mlUnusedTail = &mNodes[sz-1];
        }

        BAResult& malloc() {
            std::unique_lock<std::mutex> lk(mMutex);
            //debugInfo("req-deq");
            mCondEmpty.wait(lk, [this]{return mlUnused != InvalidPtr;});
            // transfer the node from Available list to ReqDeq list
            NodePtr cN = disconnectHead(mlUnused, mlUnusedTail);
            // connect node to mReqDeq
            append(cN, mlReserved, mlReservedTail);
            lk.unlock();
            return cN->result;
        }
        void free(uint64_t id) {
            std::lock_guard<std::mutex> lk(mMutex);
            //debugInfo("register-deq");
            // transfer node from mReqDeq to Unused
            NodePtr cN = &mNodes[id];
            // remove it from ReqEnq
            disconnectNode(cN, mlReserved, mlReservedTail);
            // connect it to Available
            append(cN, mlUnused, mlUnusedTail);
            mCondEmpty.notify_one();    
        }

        void debugInfo(const std::string& s) {
            std::cerr << std::endl << s << std::endl;
            std::cerr << "Unused:"; cerrList(mlUnused);
            std::cerr << "Rserved:";  cerrList(mlReserved);
        }

    private:

        const uint64_t mMaxSize;
        const NodePtr InvalidPtr; 

        NodePtr mlUnused;
        NodePtr mlReserved;

        NodePtr mlUnusedTail;
        NodePtr mlReservedTail;

        std::vector<Node> mNodes;
        std::mutex mMutex;
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
                if (cN->next != InvalidPtr) cN->next->prev = cN->prev;
                if (cN->prev != InvalidPtr) cN->prev->next = cN->next;
            }
        }

        NodePtr disconnectHead(NodePtr &lH, NodePtr& lT) {
            if (lH == InvalidPtr) return InvalidPtr;
            NodePtr t = lH;
            lH = lH->next;
            if (lH == InvalidPtr) lT = InvalidPtr;
            else lH->prev = InvalidPtr;
            return t;
        }
        NodePtr disconnectTail(NodePtr &lH, NodePtr& lT) {
            if (lT == InvalidPtr) return InvalidPtr;
            NodePtr t = lT;
            lT = lT->prev;
            if (lT == InvalidPtr) lH = InvalidPtr;
            else lT->next = InvalidPtr;
            return t;
        }

        void prepend(NodePtr nN, NodePtr& lH, NodePtr& lT) { 
            if (lH == InvalidPtr) {
                // empty list
                nN->prev = InvalidPtr;
                nN->next = InvalidPtr;
                lH = nN; lT = nN;
            } else {
                // there is a head already
                nN->next = lH;
                nN->prev = InvalidPtr;
                lH->prev = nN;
                lH = nN;
            }
        }  
        void append(NodePtr nN, NodePtr& lH, NodePtr& lT) { 
            if (lH == InvalidPtr) {
                // empty list so set both to this node
                lH = nN; lT = nN;
                nN->prev = InvalidPtr;
                nN->next = InvalidPtr;
            } else {
                // the list is not empty
                nN->prev = lT;
                nN->next = InvalidPtr;
                lT->next = nN;
                lT = nN;
            }
        }
};

#endif
