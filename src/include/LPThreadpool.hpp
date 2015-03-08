#ifndef __LP_THREADPOOL__
#define __LP_THREADPOOL__

#pragma once

typedef std::function<void(uint32_t, uint32_t, void*)> LPFunc;

class ISingleTaskPool {
    public:
        virtual ~ISingleTaskPool() {};
        virtual void initThreads() = 0;
        virtual void startSingleAll(LPFunc f, void *args = nullptr) = 0;
        virtual void waitSingleAll() = 0;
        virtual void destroy() = 0;
};

class IMultiTaskPool {
    public:
        virtual ~IMultiTaskPool() {};
        virtual void initThreads() = 0;
        virtual void startAll() = 0;
        virtual void waitAll() = 0;
        virtual void destroy() = 0;
        virtual void addTask(LPFunc f, void* args = nullptr) = 0;
};

#endif
