#include "QXYSquare.h"

#include <algorithm>
#include <numeric>
#include <assert.h>

#include "MathUtil.h"

using mcf::kPi;

QXYSquare::QXYSquare(int size, int Ntau, unsigned int seed)
	:
	QXYSquare(size, Ntau, 1.0, 1.0, seed)
{}

QXYSquare::QXYSquare(int size, int Ntau, float K_s, float K_t, unsigned int seed)
	:
	size(size),
	Ntau(Ntau),
	K_s(K_s),
	K_t(K_t),
	ss_size(size * size),
	st_size(size * Ntau),
	lattice({ size, size, Ntau }, { "x", "y", "tau" }),
	System(size * size * Ntau, size * size * Ntau * 3)
{
	site_fields = std::vector<double>(n_site_variables);
	plaq_fields = std::vector<double>(n_site_variables*3);

	std::mt19937 rng(seed != 0 ? seed : std::random_device{}());
	std::uniform_real_distribution<double> dst;

	std::for_each(site_fields.begin(), site_fields.end(), [&rng, &dst](double& d) {d = dst(rng); });
}

double QXYSquare::getEnergy() const
{
	double energy = 0.0;
	for (int nx = 0; nx < size; nx++)
	{
		for (int ny = 0; ny < size; ny++)
		{
			for (int nt = 0; nt < Ntau; nt++)
			{
				int siteIndex = lattice.index(nx, ny, nt);

				// No double counting
				int i_r = lattice.neighbor(siteIndex, 0, +1);
				int i_u = lattice.neighbor(siteIndex, 1, +1);
				int i_t = lattice.neighbor(siteIndex, 2, +1);

				energy += -K_s * (cos(kPi * (site_fields[i_r] - site_fields[siteIndex])) + cos(kPi * (site_fields[i_u] - site_fields[siteIndex])));
				energy += -K_t * (cos(kPi * (site_fields[i_t] - site_fields[siteIndex])));
			}
		}
	}

	return energy;
}

double QXYSquare::proposeSiteFlip(int index, double angle) const
{
	std::pair<std::vector<int>, std::vector<int>> connectedSites = getSiteConnectedCluster(index);

	double flip_energy = 0.0;
	for (const int& csite : connectedSites.first)
	{
		flip_energy += -K_s * (cos(kPi * (site_fields[index] + angle - site_fields[csite])) - cos(kPi * (site_fields[index] - site_fields[csite])));
	}
	for (const int& tsite : connectedSites.second)
	{
		double dTheta = site_fields[index] - site_fields[tsite];
		double dTheta_f = dTheta + angle;
		flip_energy += -K_t * (cos(kPi * dTheta_f) - cos(kPi * dTheta));
	}

	return flip_energy;
}

double QXYSquare::proposePlaqFlip(int index, double angle) const
{
	return 0.0;
}

void QXYSquare::UpdateSite(int index, double angle)
{
	site_fields[index] += angle;
}

void QXYSquare::UpdatePlaq(int index, double angle)
{
	plaq_fields[index] += angle;
}

std::vector<std::pair<std::vector<int>, int>> QXYSquare::getSpacialVortices() const
{
	std::vector<std::pair<std::vector<int>, int>> vortices;
	for (int n = 0; n < n_plaq_variables; n += 3)
	{
		std::vector<int> plaq_sites = getPlaqConnectedCluster(n).first;

		double vortex = 0.0;
		size_t size = plaq_sites.size();
		for (int i = 0; i < size; i++)
		{
			int site_2 = plaq_sites[(i + 1) % size];
			int site_1 = plaq_sites[i];

			double d2 = site_fields[site_2];
			double d1 = site_fields[site_1];

			vortex += mcf::mapToCircle(d2 - d1);
		}

		if (vortex >= 1.0 - 1e-5 || vortex <= -1.0 + 1e-5)
		{
			vortices.push_back(std::pair<std::vector<int>, int>(plaq_sites, mcf::sgn(vortex)));
		}
	}

	return vortices;
}

std::vector<std::pair<std::vector<int>, int>> QXYSquare::getTemporalVortices() const
{
	std::vector<std::pair<std::vector<int>, int>> vortices;
	for (int n = 0; n < 2*n_plaq_variables/3; n++)
	{
		int plaq_index = n + (n - 1) / 2; // this ignores spacial plaquettes
		std::vector<int> plaq_sites = getPlaqConnectedCluster(plaq_index).first;

		double vortex = 0.0;
		size_t size = plaq_sites.size();
		for (int i = 0; i < size; i++)
		{
			int site_2 = plaq_sites[(i + 1) % size];
			int site_1 = plaq_sites[i];

			double d2 = site_fields[site_2];
			double d1 = site_fields[site_1];

			vortex += mcf::mapToCircle(d2 - d1);
		}

		if (vortex >= 1.0 - 1e-5 || vortex <= -1.0 + 1e-5)
		{
			vortices.push_back(std::pair<std::vector<int>, int>(plaq_sites, mcf::sgn(vortex)));
		}
	}

	return vortices;
}

void QXYSquare::LogToFile(std::ofstream& outfile) const
{
	const auto vortexPairs_s = getSpacialVortices();
	const auto vortexPairs_t = getTemporalVortices();
	outfile << "(" << vortexPairs_s.size() << "," << vortexPairs_t.size() << ")";
}

std::pair<std::vector<int>, std::vector<int>> QXYSquare::getSiteConnectedCluster(int siteIndex) const
{
	int i_u = lattice.neighbor(siteIndex, 1, +1);
	int i_r = lattice.neighbor(siteIndex, 0, +1);
	int i_d = lattice.neighbor(siteIndex, 1, -1);
	int i_l = lattice.neighbor(siteIndex, 0, -1);

	int it_u = lattice.neighbor(siteIndex, 2, +1);
	int it_d = lattice.neighbor(siteIndex, 2, -1);

	return std::pair<std::vector<int>, std::vector<int>>({ i_u, i_r, i_d, i_l }, { it_u, it_d });
}

std::pair<std::vector<int>, std::vector<int>> QXYSquare::getPlaqConnectedCluster(int plaqIndex) const
{
	int plaqType = plaqIndex % 3;
	int i_a = plaqIndex / 3;

	int i_b = -1;
	int i_c = -1;
	int i_d = -1;

	switch (plaqType)
	{
	case 0: // xy
		i_b = lattice.neighbor(i_a, 1, +1);
		i_c = lattice.neighbor(i_b, 0, +1);
		i_d = lattice.neighbor(i_a, 0, +1);
		break;
	case 1: // x-tau
		i_b = lattice.neighbor(i_a, 2, +1);
		i_c = lattice.neighbor(i_b, 0, +1);
		i_d = lattice.neighbor(i_a, 0, +1);
		break;
	case 2: // y-tau
		i_b = lattice.neighbor(i_a, 2, +1);
		i_c = lattice.neighbor(i_b, 1, +1);
		i_d = lattice.neighbor(i_a, 1, +1);
		break;
	}

	return std::pair<std::vector<int>, std::vector<int>>({ i_a, i_b, i_c, i_d }, {});
}