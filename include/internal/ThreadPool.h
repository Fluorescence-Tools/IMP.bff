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
#include <chrono>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>
#ifdef __APPLE__
#include <pthread.h>
#include <sys/qos.h>
#endif

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
    std::atomic<bool> stop_{false};

    // Workers spin on the generation counter for a while before they sleep
    // on the condition variable: an evaluation issues its phases a few tens
    // of microseconds apart, and a condition-variable wake-up costs about
    // that much on its own.
    std::atomic<unsigned long> gen_atomic_{0};
    // Spin (yielding) for up to this long before sleeping: longer than the
    // gap between the phases of one evaluation and between evaluations of a
    // tight screening loop, short enough that an idle process settles.
    static constexpr double SPIN_SECONDS = 2e-3;
    static const int SPIN_ITERATIONS = 200000;

    void loop() {
        unsigned long seen = 0;
        for (;;) {
            bool got = false;
            auto t0 = std::chrono::steady_clock::now();
            for (int i = 0; i < SPIN_ITERATIONS; i++) {
                if (stop_) return;
                if (gen_atomic_.load(std::memory_order_acquire) != seen) { got = true; break; }
                std::this_thread::yield();
                if ((i & 63) == 63 &&
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() > SPIN_SECONDS) break;
            }
            if (!got) {
                std::unique_lock<std::mutex> lk(mutex_);
                start_.wait(lk, [&] { return stop_ || generation_ != seen; });
                if (stop_) return;
            }
            seen = gen_atomic_.load(std::memory_order_acquire);
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
            workers_.emplace_back([this] {
#ifdef __APPLE__
                // Ask for performance cores: the evaluation is latency-bound
                // and an efficiency core doubles a task's time.
                pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
#endif
                loop();
            });
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
            gen_atomic_.store(generation_, std::memory_order_release);
        }
        start_.notify_all();
        work();
        // the caller spins too: the workers finish within microseconds
        for (int i = 0; i < SPIN_ITERATIONS && active_.load() != 0; i++) {
            std::this_thread::yield();
        }
        if (active_.load() != 0) {
            std::unique_lock<std::mutex> lk(mutex_);
            done_.wait(lk, [&] { return active_ == 0; });
        }
        fn_ = nullptr;
    }
};

IMPBFF_END_INTERNAL_NAMESPACE

#endif /* IMPBFF_INTERNAL_THREAD_POOL_H */
