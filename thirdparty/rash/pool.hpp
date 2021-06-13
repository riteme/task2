#pragma once

#include <list>
#include <thread>
#include <future>
#include <functional>


class ThreadPool {
public:
    using TaskFn = std::function<void()>;

    ThreadPool() : ThreadPool(std::thread::hardware_concurrency()) {}
    ThreadPool(int n_workers);
    ~ThreadPool();

    ThreadPool(const ThreadPool &) = delete;
    ThreadPool(ThreadPool &&) = delete;
    auto operator=(const ThreadPool &) = delete;
    auto operator=(ThreadPool &&) = delete;

    auto run(const TaskFn &fn) -> std::future<void>;

private:
    struct Task {
        TaskFn fn;
        std::promise<void> promise;
    };

    void _worker_fn();

    std::atomic<bool> stopped;
    std::mutex mutex;
    std::condition_variable cond;
    std::vector<std::thread> workers;
    std::list<Task> tasks;
};
