#include "../include/McMachine.h"

#include <iostream>
#include <assert.h>
#include <numeric>
#include <cmath>
#include <iomanip>
#include <stdexcept>

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
	// A ground state start still draws the random configuration first, so the stream that drives
	// the sweeps is the same for both kinds of start.
	if (params.ground_state_start)
		system.setGroundState();

	logfile_name = filename;
}

void McMachine::addProvenance(std::string key, std::string value)
{
	provenance.emplace_back(std::move(key), std::move(value));
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

std::string McMachine::binsFileName(const std::string& log_name)
{
	const std::string suffix = ".txt";
	if (log_name.size() >= suffix.size()
		&& log_name.compare(log_name.size() - suffix.size(), suffix.size(), suffix) == 0)
		return log_name.substr(0, log_name.size() - suffix.size()) + "_bins.txt";
	return log_name + "_bins.txt";
}

// The '#' provenance and parameter lines, the same at the top of the log and of the bins file.
void McMachine::writeHeader(std::ofstream& out)
{
	for (const auto& entry : provenance)
		out << "# " << entry.first << '\t' << entry.second << '\n';
	out << "# seed\t" << seed << '\n';
	out << "# t_max\t" << params.t_max << "\tt_min\t" << params.t_min
		<< "\tt_fac\t" << params.t_fac << '\n';
	out << "# updates_per_sweep\t" << params.updates_per_sweep
		<< "\tinitial_therm_sweeps\t" << params.initial_therm_sweeps
		<< "\tmax_therm_sweeps\t" << params.max_therm_sweeps
		<< "\tmax_measure_sweeps\t" << params.max_measure_sweeps
		<< "\tn_measurements\t" << params.n_measurements << '\n';
	out << "# delta\t" << params.delta << "\toverrelax\t" << params.overrelax
		<< "\tupdates_per_overrelaxation\t" << params.updates_per_overrelaxation
		<< "\tground_state_start\t" << params.ground_state_start << '\n';
	out << "# variables\t" << system.numVariables() << '\n';
}

void McMachine::StartSimulation()
{
	logfile = std::ofstream(logfile_name);
	assert(logfile.is_open());
	writeHeader(logfile);

	// The column-header line stays the first non-comment line, so parsers reading with
	// comment='#' keep working. Every observable gets a value column and a 'D' error column.
	const std::vector<std::string> observable_names = system.observableNames();
	logfile << "T";
	for (const std::string& name : observable_names)
		logfile << '\t' << name << "\tD" << name;
	logfile << "\tautocorrelation\tn_sweeps\tacceptance\n";

	if (params.bin_size > 0)
	{
		binsfile = std::ofstream(binsFileName(logfile_name));
		if (!binsfile.is_open())
			throw std::runtime_error("cannot open " + binsFileName(logfile_name));
		writeHeader(binsfile);
		binsfile << "# bin_size\t" << params.bin_size << '\n';
		binsfile << "# One line per bin of consecutive measurement sweeps at one temperature, measured after\n"
			"# every sweep. <name> is the mean of the observable over the bin and <name>_sq the mean of its\n"
			"# square. bin counts from 0 at each temperature; n_sweeps is bin_size except for the last bin of a\n"
			"# temperature, which holds what is left over.\n";
		binsfile << "T\tbin\tn_sweeps";
		for (const std::string& name : observable_names)
			binsfile << '\t' << name << '\t' << name << "_sq";
		binsfile << '\n';
		binsfile << std::setprecision(17);
	}

	const int buffersize = 50;
	BufferedArray energies(buffersize);

	std::cout << "Initial Thermalization" << std::endl;

	// t_fac < 1 cools from t_max down to t_min; t_fac > 1 heats from t_min up to t_max.
	const bool heating = params.t_fac > 1.0;
	double temperature = heating ? params.t_min : params.t_max;

	Thermalize(current_nSweeps, energies, temperature);

	while (heating ? temperature < params.t_max : temperature > params.t_min)
	{
		Thermalize(current_nSweeps, energies, temperature);

		Measure(params.n_measurements, current_measurement_sweeps, temperature);

		temperature *= params.t_fac;
	}

	logfile.close();
	if (binsfile.is_open())
		binsfile.close();
}

void McMachine::writeBin(double temperature, int bin, int n_sweeps, const std::vector<double>& sum,
	const std::vector<double>& sum_sq)
{
	binsfile << temperature << '\t' << bin << '\t' << n_sweeps;
	for (size_t i = 0; i < sum.size(); i++)
		binsfile << '\t' << sum[i] / n_sweeps << '\t' << sum_sq[i] / n_sweeps;
	binsfile << '\n';
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

	// The bins: running sums of every observable and of its square over the current bin.
	const bool binning = params.bin_size > 0 && binsfile.is_open();
	const size_t n_names = binning ? system.observableNames().size() : 0;
	std::vector<double> bin_sum(n_names, 0.0), bin_sum_sq(n_names, 0.0);
	int bin = 0, in_bin = 0;

	for (int n = 0; n < n_measurements; n++)
	{
		for (int m = 0; m < n_measure_sweeps; m++)
		{
			energies.push_back(system.getEnergy());
			Sweep(params.updates_per_sweep, temperature, false);
			if (binning)
			{
				const std::vector<double> values = system.measure(temperature);
				for (size_t i = 0; i < n_names; i++)
				{
					bin_sum[i] += values[i];
					bin_sum_sq[i] += values[i] * values[i];
				}
				if (++in_bin == params.bin_size)
				{
					writeBin(temperature, bin++, in_bin, bin_sum, bin_sum_sq);
					std::fill(bin_sum.begin(), bin_sum.end(), 0.0);
					std::fill(bin_sum_sq.begin(), bin_sum_sq.end(), 0.0);
					in_bin = 0;
				}
			}
		}
		if (params.overrelax)
		{
			Overrelax(params.updates_per_overrelaxation);
		}
		samples.push_back(system.measure(temperature));
	}
	if (binning)
	{
		if (in_bin > 0)
			writeBin(temperature, bin, in_bin, bin_sum, bin_sum_sq);
		// Flushed once per temperature, like the log, so a job stopped by its time limit keeps
		// every temperature it finished.
		binsfile.flush();
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
	// The variance of a single measurement: s^2 = (1/N)Sum(x-mean)^2.
	auto error_of = [&](int column, double mean) {
		if (N < 1)
			return 0.0;
		double sum_sqr = 0.0;
		for (int n = 0; n < N; n++)
		{
			const double d = samples[n][column] - mean;
			sum_sqr += d * d;
		}
		return sum_sqr / double(N);
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