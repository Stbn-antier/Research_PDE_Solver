#pragma once
#include <thread>
#include <future>
#include <vector>
#include <queue>

class Thread_pool
{
private:
	std::queue<std::function<void()>> Tasks_;
	std::vector<std::jthread> Threads_;
	std::mutex queue_mutex_;
	std::condition_variable cv_;

public:
	Thread_pool(size_t num_threads = std::jthread::hardware_concurrency());
	void enqueue(std::function<void()> task);
};

