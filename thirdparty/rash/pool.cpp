#include "pool.hpp"


ThreadPool::ThreadPool(int n_workers)
    : stopped(false) {
    workers.reserve(n_workers);
    for (int i = 0; i < n_workers; i++) {
        auto t = std::thread([this] {
            _worker_fn();
        });

        workers.push_back(std::move(t));
    }
}

ThreadPool::~ThreadPool() {
    std::unique_lock lock(mutex);
    stopped = true;
    lock.unlock();

    cond.notify_all();

    for (auto &t : workers) {
        t.join();
    }
}

void ThreadPool::_worker_fn() {
    while (true) {
        std::unique_lock lock(mutex);
        cond.wait(lock, [this] {
            return !tasks.empty() || stopped;
        });

        if (!tasks.empty()) {
            auto task = std::move(tasks.front());
            tasks.pop_front();
            lock.unlock();

            try {
                task.fn();
                task.promise.set_value();
            } catch (std::exception_ptr e) {
                task.promise.set_exception(e);
            }
        } else
            break;
    }
}

auto ThreadPool::run(const TaskFn &fn) -> std::future<void> {
    auto task = Task{fn};
    auto future = task.promise.get_future();

    std::unique_lock lock(mutex);
    tasks.push_back(std::move(task));
    lock.unlock();

    cond.notify_one();

    return future;
}
