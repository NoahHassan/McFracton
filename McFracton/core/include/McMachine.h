#pragma once

#include <algorithm>
#include <fstream>

#include "System.h"
#include "BufferedArray.h"

class McMachine {
public:
	struct NumericalParams {
		double t_max = 100.0;
		double t_min = 0.1;
		double t_fac = 0.9;
		double delta = 0.1;
		int updates_per_sweep = 5000;
		int max_therm_sweeps = 5000;
		int max_measure_sweeps = 500;
		int n_measurements = 10;
		bool overrelax = false;
		bool log_energies = false;
		int updates_per_overrelaxation = 10000;
	};
	struct AutoCorrResult {
		double tau_int;
		std::vector<double> rho;
	};
public:
	McMachine(NumericalParams params, System& system, std::string filename);
public:
	void Sweep(int nUpdates, const float temperature);
	void Overrelax(int nUpdates);
	void StartSimulation();
private:
	void Thermalize(int maxSweeps, BufferedArray& energies, const float temperature);
	void Measure(int n_measurements, int n_measure_sweeps, const double temperature);
	AutoCorrResult Autocorrelation(const std::vector<double>& data);
private:
	NumericalParams params;
	std::mt19937 rng;
	std::uniform_int_distribution<int> site_dst;
	std::uniform_real_distribution<double> eps_dst;
	std::uniform_real_distribution<float> acc_dst;
	System& system;
	std::ofstream logfile;
	double acceptance_ratio;
	int current_nSweeps;
	int current_measurement_sweeps;
};