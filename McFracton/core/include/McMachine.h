#pragma once

#include "System.h"

class McMachine {
public:
	struct NumericalParams {
	public:
		NumericalParams() = default;
		NumericalParams(double t_max, double t_min, double eps)
			:
			t_max(t_max),
			t_min(t_min),
			eps(eps)
		{}
	public:
		double t_max;
		double t_min;
		double eps;
	};
public:
	McMachine(NumericalParams params, System& system);
public:
	void Sweep(int nSweeps, double T);
private:
	NumericalParams params;
	std::mt19937 rng;
	std::uniform_int_distribution<int> site_dst;
	std::uniform_real_distribution<double> eps_dst;
	std::uniform_real_distribution<float> acc_dst;
	System& system;
};