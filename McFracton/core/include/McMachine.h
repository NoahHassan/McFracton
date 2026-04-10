#pragma once

#include <algorithm>
#include <fstream>

#include "System.h"
#include "BufferedArray.h"

class McMachine {
public:
	struct NumericalParams {
	public:
		NumericalParams(
			double t_max = 100.0,
			double t_min = 0.1,
			double t_fac = 0.9,
			double delta = 0.1,
			int updates_per_sweep = 5000,
			int max_therm_sweeps = 5000,
			int measure_sweeps = 100
		)
			:
			t_max(t_max),
			t_min(t_min),
			t_fac(t_fac),
			delta(delta),
			max_therm_sweeps(max_therm_sweeps),
			updates_per_sweep(updates_per_sweep),
			measure_sweeps(measure_sweeps)
		{}
	public:
		double t_max;
		double t_min;
		double t_fac;
		double delta;
		int max_therm_sweeps;
		int updates_per_sweep;
		int measure_sweeps;
	};
public:
	McMachine(NumericalParams params, System& system, std::string filename);
public:
	void Sweep(int nUpdates);
	void StartSimulation();
private:
	void Thermalize(int maxSweeps, BufferedArray& energies);
	void Measure(int nSweeps);
private:
	NumericalParams params;
	std::mt19937 rng;
	std::uniform_int_distribution<int> site_dst;
	std::uniform_real_distribution<double> eps_dst;
	std::uniform_real_distribution<float> acc_dst;
	System& system;
	std::ofstream logfile;
};