#include "../include/McMachine.h"

#include <iostream>
#include <assert.h>

McMachine::McMachine(NumericalParams params, System& system, std::string filename = "log.txt")
	:
	params(params),
	system(system),
	acceptance_ratio(0.5)
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

	const int buffersize = 50;
	BufferedArray energies(buffersize);

	std::cout << "Initial Thermalization" << std::endl;

	Thermalize(params.max_therm_sweeps, energies, params.t_max);

	double temperature = params.t_max;
	while (temperature > params.t_min)
	{
		Thermalize(params.max_therm_sweeps, energies, temperature);

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

		if (n > energies.get_size())
		{
			float avg = energies.get_average();
			float dev = energies.get_deviation();

			if (dev < 1.0f * temperature)
				return;
		}
	}

	std::cout << maxSweeps << " sweeps completed" << std::endl;
}

void McMachine::Measure(int nSweeps, const double temperature)
{
	std::cout << "Measuring" << std::endl;
	logfile << temperature;
	for (int n = 0; n < nSweeps; n++)
	{
		Sweep(params.updates_per_sweep, (float)temperature);
		if (params.overrelax)
		{
			Overrelax(params.updates_per_overrelaxation);
		}
		
		logfile << '\t';
		system.LogToFile(logfile);
	}

	logfile << std::endl;
}