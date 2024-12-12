#include "thread_pool.h"

ThreadPool::ThreadPool(size_t numThreads) : stop(false), activeTasks(0) {
    pthread_mutex_init(&queueMutex, nullptr);
    pthread_cond_init(&condition, nullptr);
    pthread_mutex_init(&activeMutex, nullptr);
    pthread_cond_init(&activeCondition, nullptr);

    threads.resize(numThreads);
    for (size_t i = 0; i < numThreads; ++i) {
        pthread_create(&threads[i], nullptr, ThreadPool::worker, this);
    }
}

ThreadPool::~ThreadPool() {
    // Stop all threads
    pthread_mutex_lock(&queueMutex);
    stop = true;
    pthread_cond_broadcast(&condition);
    pthread_mutex_unlock(&queueMutex);

    // Join all threads
    for (pthread_t& thread : threads) {
        pthread_join(thread, nullptr);
    }

    pthread_mutex_destroy(&queueMutex);
    pthread_cond_destroy(&condition);
    pthread_mutex_destroy(&activeMutex);
    pthread_cond_destroy(&activeCondition);
}

void ThreadPool::enqueue(std::function<void()> task) {
    pthread_mutex_lock(&queueMutex);
    tasks.push(task);
    pthread_cond_signal(&condition);
    pthread_mutex_unlock(&queueMutex);
}

void ThreadPool::wait() {
    pthread_mutex_lock(&activeMutex);
    while (activeTasks > 0 || !tasks.empty()) {
        pthread_cond_wait(&activeCondition, &activeMutex);
    }
    pthread_mutex_unlock(&activeMutex);
}

void* ThreadPool::worker(void* arg) {
    ThreadPool* pool = static_cast<ThreadPool*>(arg);
    pool->run();
    return nullptr;
}

void ThreadPool::run() {
    while (true) {
        std::function<void()> task;

        // Acquire task
        pthread_mutex_lock(&queueMutex);
        while (tasks.empty() && !stop) {
            pthread_cond_wait(&condition, &queueMutex);
        }
        if (stop && tasks.empty()) {
            pthread_mutex_unlock(&queueMutex);
            break;
        }
        task = tasks.front();
        tasks.pop();
        pthread_mutex_unlock(&queueMutex);

        // Increment active tasks
        pthread_mutex_lock(&activeMutex);
        ++activeTasks;
        pthread_mutex_unlock(&activeMutex);

        // Execute the task
        task();

        // Decrement active tasks
        pthread_mutex_lock(&activeMutex);
        --activeTasks;
        if (activeTasks == 0 && tasks.empty()) {
            pthread_cond_signal(&activeCondition);
        }
        pthread_mutex_unlock(&activeMutex);
    }
}

