#include "../include/McMachine.h"

#include <iostream>
#include <assert.h>

McMachine::McMachine(NumericalParams params, System& system, std::string filename = "log.txt")
	:
	params(params),
	system(system)
{
	std::random_device rd;
	std::mt19937 rng(rd());
	site_dst = std::uniform_int_distribution<int>(0, system.nSites - 1);
	eps_dst = std::uniform_real_distribution<double>(-params.delta, params.delta);
	acc_dst = std::uniform_real_distribution<float>(0.0f, 1.0f);

	logfile = std::ofstream(filename);
}

void McMachine::Sweep(int nUpdates, double T)
{
	for (int n = 0; n < nUpdates; n++)
	{
		int site_index = site_dst(rng);
		double flip_angle = eps_dst(rng);
		double dE = system.proposeSiteFlip(site_index, flip_angle);

		if (dE < 0.0 || acc_dst(rng) < std::exp(-dE / T))
		{
			system.UpdateSite(site_index, flip_angle);
		}
	}
}

void McMachine::StartSimulation()
{
	assert(logfile.is_open());

	const int buffersize = 50;
	BufferedArray energies(buffersize);

	std::cout << "Initial Thermalization" << std::endl;

	Thermalize(params.t_max, params.max_therm_sweeps, energies);

	double temperature = params.t_max;
	while (temperature > params.t_min)
	{
		Thermalize(temperature, params.max_therm_sweeps, energies);

		Measure(params.measure_sweeps, temperature);

		temperature *= params.t_fac;
	}

	logfile.close();
}

void McMachine::Thermalize(double T, int maxSweeps, BufferedArray& energies)
{
	std::cout << "Thermalizing at T = " << T << std::endl;
	for (int n = 0; n < maxSweeps; n++)
	{
		Sweep(params.updates_per_sweep, T);
		energies.Push((float)system.getEnergy());

		if (n > energies.get_size())
		{
			float avg = energies.get_average();
			float dev = energies.get_deviation();

			if (dev < 1.0f * T)
				return;
		}
	}
}

void McMachine::Measure(int nSweeps, double T)
{
	logfile << T;

	std::cout << "Measuring" << std::endl;
	for (int n = 0; n < nSweeps; n++)
	{
		Sweep(params.updates_per_sweep, T);
		
		logfile << '\t';
		system.LogToFile(logfile);
	}

	logfile << std::endl;
}