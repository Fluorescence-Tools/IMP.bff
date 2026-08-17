/**
 *  \file IMP/bff/internal/ThreadPool.h
 *  \brief A minimal persistent thread pool for the AV evaluation phases.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_INTERNAL_THREAD_POOL_H
#define IMPBFF_INTERNAL_THREAD_POOL_H

#include <IMP/bff/bff_config.h>

#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

IMPBFF_BEGIN_INTERNAL_NAMESPACE

//! Persistent workers that run `fn(i)` for i in [0, n) with dynamic
//! scheduling; run() blocks until every task finished. The calling thread
//! works too. Not reentrant.
class ThreadPool {
    std::vector<std::thread> workers_;
    std::mutex mutex_;
    std::condition_variable start_, done_;
    unsigned long generation_ = 0;
    std::atomic<size_t> next_{0};
    size_t n_tasks_ = 0;
    std::atomic<int> active_{0};
    const std::function<void(size_t)> *fn_ = nullptr;
    bool stop_ = false;

    void loop() {
        unsigned long seen = 0;
        for (;;) {
            {
                std::unique_lock<std::mutex> lk(mutex_);
                start_.wait(lk, [&] { return stop_ || generation_ != seen; });
                if (stop_) return;
                seen = generation_;
            }
            work();
            if (--active_ == 0) {
                std::lock_guard<std::mutex> lk(mutex_);
                done_.notify_all();
            }
        }
    }
    void work() {
        for (;;) {
            size_t i = next_++;
            if (i >= n_tasks_) return;
            (*fn_)(i);
        }
    }

public:
    explicit ThreadPool(int n_workers) {
        for (int i = 1; i < n_workers; i++) {   // the caller is worker 0
            workers_.emplace_back([this] { loop(); });
        }
    }
    ~ThreadPool() {
        {
            std::lock_guard<std::mutex> lk(mutex_);
            stop_ = true;
        }
        start_.notify_all();
        for (auto &w : workers_) w.join();
    }
    int size() const { return (int) workers_.size() + 1; }

    void run(size_t n, const std::function<void(size_t)> &fn) {
        if (n == 0) return;
        fn_ = &fn;
        n_tasks_ = n;
        next_ = 0;
        active_ = (int) workers_.size();
        {
            std::lock_guard<std::mutex> lk(mutex_);
            generation_++;
        }
        start_.notify_all();
        work();
        std::unique_lock<std::mutex> lk(mutex_);
        done_.wait(lk, [&] { return active_ == 0; });
        fn_ = nullptr;
    }
};

IMPBFF_END_INTERNAL_NAMESPACE

#endif /* IMPBFF_INTERNAL_THREAD_POOL_H */
