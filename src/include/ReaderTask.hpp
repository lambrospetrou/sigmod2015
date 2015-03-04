#ifndef __LP_READER_TASK__
#define __LP_READER_TASK__

#pragma once

#include <iostream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "LPUtils.hpp"
#include "ReferenceTypes.hpp"

namespace lp {

    struct ReceivedMessage {
        MessageHead head;
        std::vector<char> data;
    };
    
    class ReaderIO {
        public:
            static ReceivedMessage* readInput(FILE *instream);
            static ReceivedMessage* readInput(std::istream& instream);
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
};

#endif
