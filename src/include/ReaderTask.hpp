#ifndef __LP_READER_TASK__
#define __LP_READER_TASK__

#pragma once

#include <iostream>
#include <ios>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <thread>
#include <chrono>

#include "LPUtils.hpp"
#include "ReferenceTypes.hpp"


#include "circularfifo_memory_relaxed_aquire_release.hpp"
#include "concurrentqueue.h"
#include "readerwriterqueue.h"

namespace lp {

    struct ReceivedMessage {
        MessageHead head;
        std::vector<char> data;
    };

    //typedef CircularFifo<ReceivedMessage*, 5000> LPMsgQ;
    //typedef moodycamel::ConcurrentQueue<ReceivedMessage*> LPMsgQ;
    typedef moodycamel::ReaderWriterQueue<ReceivedMessage*> LPMsgQ;

    class ReaderIO {
        public:
            ReaderIO() {
                std::ios_base::sync_with_stdio(false);
                setvbuf ( stdin , NULL , _IOFBF , 1<<22 );
            }
            virtual ~ReaderIO() {};
            static ReceivedMessage* readInput(FILE *instream);
            static ReceivedMessage* readInput(std::istream& instream);
            virtual ReceivedMessage* nextMsg() = 0;
    };

    ReceivedMessage* ReaderIO::readInput(FILE *instream) { 
        // request place from the message queue - it blocks if full
        ReceivedMessage *msg = new ReceivedMessage();
        auto& head = msg->head;
        auto& msgData = msg->data;
        // read the head of the message - type and len
        // Read the message body and cast it to the desired type
        size_t rd = fread(reinterpret_cast<char*>(&head), sizeof(head), 1, instream);
        if (unlikely(rd < 1)) { std::cerr << "read error" << std::endl; abort(); } // crude error handling, should never happen

        if (unlikely(head.type == MessageHead::Done)) {
            // exit the loop since the reader has finished its job
            return msg;
        }

        // read the actual message content
        if (unlikely(head.messageLen > msgData.size())) {
            msgData.reserve(head.messageLen);
            msgData.resize(head.messageLen);
        }

        // read the actual message content
        rd = fread(reinterpret_cast<char*>(msgData.data()), 1, head.messageLen, instream);
        // crude error handling, should never happen
        if (unlikely(rd < head.messageLen)) { std::cerr << "read error" << std::endl; abort(); }                 
        return msg;
    }


    ReceivedMessage* ReaderIO::readInput(std::istream& instream) { 
        // request place from the message queue - it blocks if full
        ReceivedMessage *msg = new ReceivedMessage();
        auto& head = msg->head;
        auto& msgData = msg->data;
        // read the head of the message - type and len
        // Read the message body and cast it to the desired type
        instream.read(reinterpret_cast<char*>(&head),sizeof(head));
        if (unlikely(!instream)) { std::cerr << "read error" << std::endl; abort(); } // crude error handling, should never happen
        if (unlikely(head.type == MessageHead::Done)) {
            // exit the loop since the reader has finished its job
            return msg;
        }
        // read the actual message content
        if (unlikely(head.messageLen > msgData.size())) {
            msgData.reserve(head.messageLen);
            msgData.resize(head.messageLen);
        }
        // read the actual message content
        instream.read(msgData.data(), head.messageLen);
        if (unlikely(!instream)) { std::cerr << "read error" << std::endl; abort(); }                 
        return msg;
    }


    class CReaderIO : public ReaderIO {
        private:
            FILE *mInstream;
        public:
            CReaderIO(FILE *instream) : ReaderIO(), mInstream(instream) {}
            ~CReaderIO() {}
            inline ReceivedMessage* __attribute__((always_inline)) nextMsg() override {
                return ReaderIO::readInput(mInstream);   
            }
    };
    class CppReaderIO : public ReaderIO {
        private:
            std::istream& mInstream;
        public:
            CppReaderIO(std::istream& instream) : ReaderIO(), mInstream(instream) {}
            ~CppReaderIO(){}
            virtual inline ReceivedMessage* __attribute__((always_inline)) nextMsg() override {
                return ReaderIO::readInput(mInstream);   
            }
    };

    class CppReaderIOTest : public ReaderIO {
        private:
            std::istream& mInstream;
            int32_t mLen;
            int32_t mEnd;
        public:
            CppReaderIOTest(std::istream& instream) : ReaderIO(), mInstream(instream) {
                mInstream.read((char*)&mLen, sizeof(int32_t));
                if (mInstream.eof()) { std::cerr << "EOF while reading file!!!" << std::endl; abort();};
                mEnd = (int32_t)mInstream.tellg() + mLen;
                //std::cerr << "len:" << mLen << std::endl;
            }
            ~CppReaderIOTest() {}
            virtual inline ReceivedMessage* __attribute__((always_inline)) nextMsg() override {
                if (mLen > 0 && (int32_t)mInstream.tellg()==mEnd) {
                    //Skip result reference data
                    mInstream.read((char*)&mLen, sizeof(int32_t));
                    mInstream.seekg(mLen, std::ios_base::cur);
                    mInstream.read((char*)&mLen, sizeof(int32_t));
                    if (mInstream.eof()) { std::cerr << "EOF while reading file!!!" << std::endl; return nullptr; }
                    mEnd = (int32_t)mInstream.tellg() + mLen;
                    //std::cerr << "len:" << mLen << std::endl;
                }
                return ReaderIO::readInput(mInstream);
            }
    };

    class AsyncReaderIO : public ReaderIO {
        private: 
            ReaderIO *mReader;
            LPMsgQ *mMsgQ;
            std::thread *mThread;
        public:
            AsyncReaderIO(ReaderIO *reader) : mReader(reader) {
                mMsgQ = new LPMsgQ();
                mThread = new std::thread([this] { ReaderTask(mMsgQ, mReader); });
            }
            ~AsyncReaderIO() {  
                mThread->join();
                delete mThread;
                delete mReader;
                delete mMsgQ;
            }
            
            // this will be called by an outsider - Single Reader ONLY
            inline ReceivedMessage* __attribute__((always_inline)) nextMsg() override {
                ReceivedMessage *msg; 
                while (!mMsgQ->try_dequeue(msg)) {/*cerr<<"m ";*/ lp::lp_spin_sleep(std::chrono::microseconds(0));}
                return msg;   
            }
        
        protected:
            // this is the Single Producer on the other side of the queue
            void ReaderTask(LPMsgQ *msgQ, ReaderIO *msgReader) {
#ifdef LPDEBUG
                //auto start = LPTimer.getChrono();
#endif
                while (true) {
                    //ReceivedMessage *msg = ReaderIO::readInput(stdin);
                    ReceivedMessage *msg = msgReader->nextMsg();
                    if (lp::unlikely(msg->head.type == MessageHead::Done)) {
                        // exit the loop since the reader has finished its job
                        //while (!msgQ.push(msg)) { lp_spin_sleep(); }
                        while (!msgQ->enqueue(msg)) { lp::lp_spin_sleep(); }
#ifdef LPDEBUG
                        //LPTimer.readingTotal += LPTimer.getChrono(start);
#endif
                        return;
                    }
                    //while (!msgQ.push(msg)) { /*cerr << "r" << std::endl;*/ lp_spin_sleep(std::chrono::microseconds(0)); }
                    while (!msgQ->enqueue(msg)) { /*cerr << "r" << std::endl;*/ lp::lp_spin_sleep(std::chrono::microseconds(0)); }
                }
                return;
            }
    };

    // Factory class that creates the proper ReaderIO class based on the arguments
    class ReaderIOFactory {
        public:
            static ReaderIO* create(FILE *mCin, bool simulateTestdriver);
            static ReaderIO* create(std::istream& mCin, bool simulateTestdriver);
            
            static ReaderIO* createAsync(FILE *mCin, bool simulateTestdriver);
            static ReaderIO* createAsync(std::istream& mCin, bool simulateTestdriver);
    };

    ReaderIO* ReaderIOFactory::create(FILE *mCin, bool simulateTestdriver = false) {
        (void)simulateTestdriver;
        return new CReaderIO(mCin);
    }
    ReaderIO* ReaderIOFactory::create(std::istream& mCin, bool simulateTestdriver = false) {
        if (simulateTestdriver) return new CppReaderIOTest(mCin);
        return new CppReaderIO(mCin);
    }
    ReaderIO* ReaderIOFactory::createAsync(FILE *mCin, bool simulateTestdriver = false) {
        return new AsyncReaderIO(create(mCin, simulateTestdriver));
    }
    ReaderIO* ReaderIOFactory::createAsync(std::istream& mCin, bool simulateTestdriver = false) {
        return new AsyncReaderIO(create(mCin, simulateTestdriver));
    }

};

#endif
