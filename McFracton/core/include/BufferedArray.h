#pragma once

#include <ranges>

class BufferedArray {
public:
	BufferedArray(int size)
		:
		max_size(size)
	{
		values.resize(max_size, 0.0f);
	}
public:
	void Push(float v) {
		values[offset] = v;
		offset = (offset + 1) % max_size;
	}
	float get_min() const {
		return *std::min_element(values.begin(), values.end());
	}
	float get_max() const {
		return *std::max_element(values.begin(), values.end());
	}
	float get_average() const {
		float sum = 0.0f;
		std::for_each(values.begin(), values.end(),
			[&sum](float v) {
				sum += v;
			}
		);
		return sum / (float)max_size;
	}
	float get_deviation() const {
		float avg = get_average();
		float sigma = 0.0f;
		std::for_each(values.begin(), values.end(),
			[&sigma, &avg](float v) {
				sigma += (v - avg) * (v - avg);
			}
		);
		return sigma / float(max_size * max_size);
	}
	int get_size() const {
		return max_size;
	}
	int get_offset() const {
		return offset;
	}
	const std::vector<float>& get_data() const {
		return values;
	}
private:
	std::vector<float> values;
	int max_size;
	int offset = 0;
};