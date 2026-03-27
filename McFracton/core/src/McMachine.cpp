#include "../include/McMachine.h"

McMachine::McMachine(NumericalParams params, System& system)
	:
	params(params),
	system(system)
{
	std::random_device rd;
	std::mt19937 rng(rd());
	site_dst = std::uniform_int_distribution<int>(0, system.nSites - 1);
	eps_dst = std::uniform_real_distribution<double>(-params.eps, params.eps);
	acc_dst = std::uniform_real_distribution<float>(0.0f, 1.0f);
}

void McMachine::Sweep(int nSweeps, double T)
{
	for (int n = 0; n < nSweeps; n++)
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