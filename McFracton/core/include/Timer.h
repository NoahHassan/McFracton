#pragma once

#include <chrono>

class Timer
{
public:
	Timer()
	{
		t0 = std::chrono::steady_clock::now();
	}
	// Seconds since construction or the last reset().
	float elapsed() const
	{
		const std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
		return std::chrono::duration<float>(t1 - t0).count();
	}
	void reset()
	{
		t0 = std::chrono::steady_clock::now();
	}
private:
	std::chrono::steady_clock::time_point t0;
};