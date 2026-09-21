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
	rng = std::mt19937(seed != 0 ? seed : std::random_device{}());
	site_dst = std::uniform_int_distribution<int>(0, system.n_site_variables - 1);
	eps_dst = std::uniform_real_distribution<double>(-1.0, 1.0);
	acc_dst = std::uniform_real_distribution<double>(0.0, 1.0);

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
		double dE = system.proposeSiteFlip(site_index, flip_angle);

		double fac = std::exp(-dE / temperature);

		if (dE < 0.0 || acc_dst(rng) < fac)
		{
			system.UpdateSite(site_index, flip_angle);
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
		system.OverrelaxSite(site_index);
	}
}

void McMachine::StartSimulation()
{
	logfile = std::ofstream(logfile_name);
	assert(logfile.is_open());
	logfile << "T\tEnergy\tDE\tHelicity Modulus\tDHM\tdefects_a\tDna\tdefects_b\tDnb\tPolyakov Loop\tDPL\tautocorrelation\tn_sweeps\tacceptance\n";

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

	std::vector<System::Observables> observables_T;
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
		observables_T.push_back(system.Measure(temperature));
	}

	// Compute observables
	const int N = (int)observables_T.size();

	auto mean_of = [&](auto getter) {
		double sum = 0.0;
		for (int n = 0; n < N; n++)
			sum += (double)getter(observables_T[n]);
		return sum / N;
		};
	auto error_of = [&](auto getter, double mean) {
		if (N < 2)
			return 0.0;
		double sum_sqr = 0.0;
		for (int n = 0; n < N; n++)
		{
			const double d = (double)getter(observables_T[n]) - mean;
			sum_sqr += d * d;
		}
		return std::sqrt(sum_sqr / (double(N) * double(N - 1)));
		};

	const double mean_energy = mean_of([](const System::Observables& o) { return o.energy; });
	const double mean_helicity = mean_of([](const System::Observables& o) { return o.helicity_modulus; });
	const double mean_defects_a = mean_of([](const System::Observables& o) { return o.n_defects_a; });
	const double mean_defects_b = mean_of([](const System::Observables& o) { return o.n_defects_b; });
	const double mean_polyakov = mean_of([](const System::Observables& o) { return o.polyakov_loop; });

	// Log observables
	logfile << temperature;
	logfile << '\t' << mean_energy << '\t' << error_of([](const System::Observables& o) { return o.energy; }, mean_energy);
	logfile << '\t' << mean_helicity << '\t' << error_of([](const System::Observables& o) { return o.helicity_modulus; }, mean_helicity);
	logfile << '\t' << mean_defects_a << '\t' << error_of([](const System::Observables& o) { return o.n_defects_a; }, mean_defects_a);
	logfile << '\t' << mean_defects_b << '\t' << error_of([](const System::Observables& o) { return o.n_defects_b; }, mean_defects_b);
	logfile << '\t' << mean_polyakov << '\t' << error_of([](const System::Observables& o) { return o.polyakov_loop; }, mean_polyakov);

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