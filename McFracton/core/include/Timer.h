#pragma once

#include <chrono>

class Timer
{
public:
	Timer()
	{
		t0 = std::chrono::steady_clock::now();
	}
	float elapsed() const
	{
		std::chrono::steady_clock::time_point t1;
		t1 = std::chrono::steady_clock::now();
		auto dt = t1 - t0;
		return (float)dt.count();
	}
	void reset()
	{
		t0 = std::chrono::steady_clock::now();
	}
private:
	std::chrono::steady_clock::time_point t0;
};