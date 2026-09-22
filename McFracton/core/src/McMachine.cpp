#include "../include/McMachine.h"

#include <iostream>
#include <assert.h>
#include <numeric>
#include <cmath>

McMachine::McMachine(NumericalParams params, System& system, std::string filename, unsigned int seed)
	:
	params(params),
	system(system),
	acceptance_ratio(0.5),
	current_nSweeps(params.initial_therm_sweeps),
	current_measurement_sweeps(params.initial_therm_sweeps)
{
	// seed == 0 means "pick a fresh, unpredictable seed"; any other value is reproducible.
	this->seed = seed != 0 ? seed : std::random_device{}();
	rng = std::mt19937(this->seed);
	site_dst = std::uniform_int_distribution<int>(0, system.numVariables() - 1);
	eps_dst = std::uniform_real_distribution<double>(-1.0, 1.0);
	acc_dst = std::uniform_real_distribution<double>(0.0, 1.0);

	// The initial configuration is the first thing drawn from the stream, before any sweep.
	system.randomize(rng);

	logfile_name = filename;
}

void McMachine::Sweep(int nUpdates, const double temperature, bool adapt_step)
{
	if (adapt_step)
		params.delta = std::max(1e-10, std::min(params.delta / (2.0 * (1.0 - acceptance_ratio)), 1.0));
	int n_accept = 0;
	for (int n = 0; n < nUpdates; n++)
	{
		int site_index = site_dst(rng);
		double flip_angle = eps_dst(rng) * params.delta;
		double dE = system.proposeUpdate(site_index, flip_angle);

		double fac = std::exp(-dE / temperature);

		if (dE < 0.0 || acc_dst(rng) < fac)
		{
			system.applyUpdate(site_index, flip_angle);
			n_accept++;
		}
	}
	acceptance_ratio = double(n_accept) / double(nUpdates);
	//std::cout << acceptance_ratio << "\t" << params.delta << std::endl;
}

void McMachine::Overrelax(int nUpdates)
{
	for (int n = 0; n < nUpdates; n++)
	{
		int site_index = site_dst(rng);
		system.overrelax(site_index, rng);
	}
}

void McMachine::StartSimulation()
{
	logfile = std::ofstream(logfile_name);
	assert(logfile.is_open());
	logfile << "# seed\t" << seed << '\n';
	logfile << "# t_max\t" << params.t_max << "\tt_min\t" << params.t_min
		<< "\tt_fac\t" << params.t_fac << '\n';
	logfile << "# updates_per_sweep\t" << params.updates_per_sweep
		<< "\tinitial_therm_sweeps\t" << params.initial_therm_sweeps
		<< "\tmax_therm_sweeps\t" << params.max_therm_sweeps
		<< "\tmax_measure_sweeps\t" << params.max_measure_sweeps
		<< "\tn_measurements\t" << params.n_measurements << '\n';
	logfile << "# delta\t" << params.delta << "\toverrelax\t" << params.overrelax
		<< "\tupdates_per_overrelaxation\t" << params.updates_per_overrelaxation << '\n';
	logfile << "# variables\t" << system.numVariables() << '\n';

	// The column-header line stays the first non-comment line, so parsers reading with
	// comment='#' keep working. Every observable gets a value column and a 'D' error column.
	const std::vector<std::string> observable_names = system.observableNames();
	logfile << "T";
	for (const std::string& name : observable_names)
		logfile << '\t' << name << "\tD" << name;
	logfile << "\tautocorrelation\tn_sweeps\tacceptance\n";

	const int buffersize = 50;
	BufferedArray energies(buffersize);

	std::cout << "Initial Thermalization" << std::endl;

	Thermalize(current_nSweeps, energies, params.t_max);

	double temperature = params.t_max;
	while (temperature > params.t_min)
	{
		Thermalize(current_nSweeps, energies, temperature);

		Measure(params.n_measurements, current_measurement_sweeps, temperature);

		temperature *= params.t_fac;
	}

	logfile.close();
}

void McMachine::Thermalize(int maxSweeps, BufferedArray& energies, const double temperature)
{
	std::cout << "Thermalizing at T = " << temperature << std::endl;
	for (int n = 0; n < maxSweeps; n++)
	{
		Sweep(params.updates_per_sweep, temperature);
		if (params.overrelax)
		{
			Overrelax(params.updates_per_overrelaxation);
		}
		if(params.log_energies)
			energies.Push((float)system.getEnergy());
	}

	std::cout << maxSweeps << " sweeps completed" << std::endl;
}

void McMachine::Measure(int n_measurements, int n_measure_sweeps, const double temperature)
{
	std::cout << "Measuring" << std::endl;

	// One row per measurement, each row holding the system's observables in declaration order.
	std::vector<std::vector<double>> samples;
	std::vector<double> energies;
	for (int n = 0; n < n_measurements; n++)
	{
		for (int m = 0; m < n_measure_sweeps; m++)
		{
			energies.push_back(system.getEnergy());
			Sweep(params.updates_per_sweep, temperature, false);
		}
		if (params.overrelax)
		{
			Overrelax(params.updates_per_overrelaxation);
		}
		samples.push_back(system.measure(temperature));
	}

	// Compute observables
	const int N = (int)samples.size();
	const int n_observables = N > 0 ? (int)samples[0].size() : 0;

	auto mean_of = [&](int column) {
		double sum = 0.0;
		for (int n = 0; n < N; n++)
			sum += samples[n][column];
		return sum / N;
		};
	// The error bar of the mean: s/sqrt(N), with s^2 using the N-1 convention.
	auto error_of = [&](int column, double mean) {
		if (N < 2)
			return 0.0;
		double sum_sqr = 0.0;
		for (int n = 0; n < N; n++)
		{
			const double d = samples[n][column] - mean;
			sum_sqr += d * d;
		}
		return std::sqrt(sum_sqr / (double(N) * double(N - 1)));
		};

	// Log observables
	logfile << temperature;
	for (int column = 0; column < n_observables; column++)
	{
		const double mean = mean_of(column);
		logfile << '\t' << mean << '\t' << error_of(column, mean);
	}

	// compute autocorrelation
	auto ac = Autocorrelation(energies);

	int required = static_cast<int>(std::max(10, (int)std::ceil(2.0 * ac.tau_int * 10)));

	current_nSweeps = std::min(params.max_therm_sweeps, required);
	current_measurement_sweeps = std::min(params.max_measure_sweeps, required);
	std::cout << "required sweeps: " << required << ", setting nSweeps = " << current_nSweeps << std::endl;

	logfile << '\t' << required << '\t' << current_nSweeps << '\t' << acceptance_ratio;
	logfile << std::endl;
}

McMachine::AutoCorrResult McMachine::Autocorrelation(const std::vector<double>& data)
{
	const size_t N = data.size();

	double mean = std::accumulate(data.begin(), data.end(), 0.0) / N;

	double var = 0.0;
	for (double x : data)
		var += (x - mean) * (x - mean);
	var /= N;

	std::vector<double> rho(N);

	rho[0] = 1.0;

	for (int t = 1; t < N; ++t) {
		double c = 0.0;

		for (int i = 0; i < N - t; ++i)
			c += (data[i] - mean) * (data[i + t] - mean);

		c /= (N - t);

		rho[t] = c / var;
	}

	// Self-consistent windowing
	double tau = 0.5;

	for (int t = 1; t < N; ++t) {
		tau += rho[t];

		if (t > 5.0 * tau)
			break;
	}
	
	return { tau, rho };
}