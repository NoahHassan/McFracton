#pragma once

#include <algorithm>
#include <fstream>
#include <random>
#include <string>
#include <utility>
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
		bool ground_state_start = false;
		// > 0: during the measurement phase, measure after every sweep as well and write the means of
		// every observable and of its square over blocks of bin_size sweeps to a second file,
		// <log>_bins.txt (see binsFileName). The chain and the main log are the same either way:
		// measuring draws no random numbers.
		int bin_size = 0;
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
	// Extra '#' provenance lines, written above the parameter lines in that order. McMachine cannot
	// know the system name, the constructor parameters or the git hash, so the caller supplies them.
	void addProvenance(std::string key, std::string value);
	void Sweep(int nUpdates, const double temperature, bool adapt_step = true);
	void Overrelax(int nUpdates);
	void StartSimulation();
	// The bins file that goes with a log: name.txt -> name_bins.txt, anything else gets _bins.txt appended.
	static std::string binsFileName(const std::string& log_name);
private:
	void writeHeader(std::ofstream& out);
	void writeBin(double temperature, int bin, int n_sweeps, const std::vector<double>& sum, const std::vector<double>& sum_sq);
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
	std::vector<std::pair<std::string, std::string>> provenance;
	std::string logfile_name;
	std::ofstream logfile;
	std::ofstream binsfile;
	double acceptance_ratio;
	int current_nSweeps;
	int current_measurement_sweeps;
};