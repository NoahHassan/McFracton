#pragma once

#include <algorithm>
#include <fstream>
#include <random>
#include <string>
#include <vector>

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
		int initial_therm_sweeps = 100;
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
	// Owns the one generator of the run. The constructor seeds it and immediately draws the
	// system's initial configuration from it, so a single seed reproduces the whole run.
	McMachine(NumericalParams params, System& system, std::string filename = "log.txt", unsigned int seed = 0);
public:
	void Sweep(int nUpdates, const double temperature, bool adapt_step = true);
	void Overrelax(int nUpdates);
	void StartSimulation();
private:
	void Thermalize(int maxSweeps, BufferedArray& energies, const double temperature);
	void Measure(int n_measurements, int n_measure_sweeps, const double temperature);
	AutoCorrResult Autocorrelation(const std::vector<double>& data);
private:
	NumericalParams params;
	// The seed actually used, so a log file records what would reproduce the run.
	unsigned int seed;
	std::mt19937 rng;
	std::uniform_int_distribution<int> site_dst;
	std::uniform_real_distribution<double> eps_dst;
	std::uniform_real_distribution<double> acc_dst;
	System& system;
	std::string logfile_name;
	std::ofstream logfile;
	double acceptance_ratio;
	int current_nSweeps;
	int current_measurement_sweeps;
};