#include "Thread_pool.h"

Thread_pool::Thread_pool(size_t num_threads)
{
    for (size_t i = 0; i < num_threads; i++)
    {
        Threads_.emplace_back(
            [this]()
            {
                while (true)
                {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(queue_mutex_);
                        cv_.wait(lock, [this]()
                            { return !Tasks_.empty(); });
                        if (Tasks_.empty())
                            return;
                        task = std::move(Tasks_.front());
                        Tasks_.pop();
                        // unlock before executing so that other tasks can enqueue task
                    }
                    task();
                }
            });
    }

}

void Thread_pool::enqueue(std::function<void()> task)
{
	std::unique_lock<std::mutex> lock(queue_mutex_);
	Tasks_.emplace(std::move(task));
	cv_.notify_one();
}

