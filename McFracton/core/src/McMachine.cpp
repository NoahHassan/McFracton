#include "../include/McMachine.h"

#include <iostream>
#include <assert.h>
#include <numeric>

McMachine::McMachine(NumericalParams params, System& system, std::string filename = "log.txt")
	:
	params(params),
	system(system),
	acceptance_ratio(0.5),
	current_nSweeps(params.max_therm_sweeps / 10)
{
	std::random_device rd;
	rng = std::mt19937(rd());
	site_dst = std::uniform_int_distribution<int>(0, system.n_site_variables - 1);
	eps_dst = std::uniform_real_distribution<double>(-1.0, 1.0);
	acc_dst = std::uniform_real_distribution<float>(0.0f, 1.0f);

	logfile = std::ofstream(filename);
}

void McMachine::Sweep(int nUpdates, const float temperature)
{
	params.delta = std::max(1e-4, std::min(params.delta / (2.0 * (1.0 - acceptance_ratio)), 1.0));
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
	assert(logfile.is_open());
	logfile << "T\tEnergy\tDE\tHelicity Modulus\tDHM\tdefects_a\tDna\tdefects_b\tDnb\tPolyakovloop\tDPL\n";

	const int buffersize = 50;
	BufferedArray energies(buffersize);

	std::cout << "Initial Thermalization" << std::endl;

	Thermalize(params.max_therm_sweeps, energies, params.t_max);

	double temperature = params.t_max;
	while (temperature > params.t_min)
	{
		Thermalize(current_nSweeps, energies, temperature);

		Measure(params.measure_sweeps, temperature);

		temperature *= params.t_fac;
	}

	logfile.close();
}

void McMachine::Thermalize(int maxSweeps, BufferedArray& energies, const float temperature)
{
	std::cout << "Thermalizing at T = " << temperature << std::endl;
	for (int n = 0; n < maxSweeps; n++)
	{
		Sweep(params.updates_per_sweep, temperature);
		if (params.overrelax)
		{
			Overrelax(params.updates_per_overrelaxation);
		}
		energies.Push((float)system.getEnergy());

		//if (n > energies.get_size())
		//{
		//	float avg = energies.get_average();
		//	float dev = energies.get_deviation();

		//	if (dev < 1.0f * temperature)
		//		return;
		//}
	}

	std::cout << maxSweeps << " sweeps completed" << std::endl;
}

void McMachine::Measure(int nSweeps, const double temperature)
{
	std::cout << "Measuring" << std::endl;

	std::vector<System::Observables> observables_T;
	for (int n = 0; n < nSweeps; n++)
	{
		Sweep(params.updates_per_sweep, (float)temperature);
		if (params.overrelax)
		{
			//Overrelax(params.updates_per_overrelaxation);
		}

		observables_T.push_back(system.Measure());
	}

	// Compute observables
	System::Observables means;
	System::Observables s_sqr;
	int N = (int)observables_T.size();
	std::vector<double> energies(N);
	for (int n = 0; n < N; n++)
	{
		energies[n] = observables_T[n].energy;

		means.energy += observables_T[n].energy / N;
		means.helicity_modulus += observables_T[n].helicity_modulus / N;
		means.n_defects_a += observables_T[n].n_defects_a;
		means.n_defects_b += observables_T[n].n_defects_b;
		means.polyakov_loop += observables_T[n].polyakov_loop / N;
	}
	means.n_defects_a /= N;
	means.n_defects_b /= N;

	for (int n = 0; n < N; n++)
	{
		s_sqr.energy += std::powf((observables_T[n].energy - means.energy), 2.0f) / N;
		s_sqr.helicity_modulus += std::powf((observables_T[n].helicity_modulus - means.helicity_modulus), 2.0f) / N;
		s_sqr.n_defects_a += std::powf(((float)observables_T[n].n_defects_a - means.n_defects_a), 2.0f);
		s_sqr.n_defects_b += std::powf(((float)observables_T[n].n_defects_b - means.n_defects_b), 2.0f);
		s_sqr.polyakov_loop += std::powf(((float)observables_T[n].polyakov_loop - means.polyakov_loop), 2.0f) / N;
	}
	s_sqr.n_defects_a /= N;
	s_sqr.n_defects_b /= N;

	// Log observables
	logfile << temperature;
	logfile << '\t' << means.energy << '\t' << s_sqr.energy;
	logfile << '\t' << means.helicity_modulus << '\t' << s_sqr.helicity_modulus;
	logfile << '\t' << means.n_defects_a << '\t' << s_sqr.n_defects_a;
	logfile << '\t' << means.n_defects_b << '\t' << s_sqr.n_defects_b;
	logfile << '\t' << means.polyakov_loop << '\t' << s_sqr.polyakov_loop;

	logfile << std::endl;

	// compute autocorrelation
	auto ac = Autocorrelation(energies);

	int required = static_cast<int>(std::ceil(2.0 * ac.tau_int * 200));

	current_nSweeps = std::min(params.max_therm_sweeps, required);
	std::cout << "required sweeps: " << required << ", setting nSweeps = " << current_nSweeps << std::endl;
}

McMachine::AutoCorrResult McMachine::Autocorrelation(const std::vector<double>& data)
{
	const int N = data.size();

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