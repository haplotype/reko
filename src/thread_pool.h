#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <pthread.h>
#include <queue>
#include <functional>
#include <vector>

class ThreadPool {
public:
    ThreadPool(size_t numThreads);
    ~ThreadPool();

    void enqueue(std::function<void()> task);
    void wait();

private:
    static void* worker(void* arg);
    void run();

    std::vector<pthread_t> threads;
    std::queue<std::function<void()>> tasks;
    pthread_mutex_t queueMutex;
    pthread_cond_t condition;
    bool stop;
    size_t activeTasks;
    pthread_mutex_t activeMutex;
    pthread_cond_t activeCondition;
};

#endif // THREADPOOL_H

