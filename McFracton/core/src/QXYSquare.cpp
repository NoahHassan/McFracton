#include "QXYSquare.h"

#include <algorithm>
#include <numeric>
#include <random>
#include <assert.h>

#include "MathUtil.h"

using mcf::kPi;

QXYSquare::QXYSquare(int size, int Ntau, unsigned int seed)
	:
	QXYSquare(size, Ntau, 1.0, 1.0, seed)
{}

QXYSquare::QXYSquare(int size, int Ntau, float K_s, float K_t, unsigned int seed)
	:
	System(mcf::PeriodicLattice({ size, size, Ntau }, { "x", "y", "tau" }), size * size * Ntau),
	size(size),
	Ntau(Ntau),
	ss_size(size * size),
	st_size(size * Ntau),
	nPlaqs(size * size * Ntau * 3),
	K_s(K_s),
	K_t(K_t)
{
	std::mt19937 rng(seed != 0 ? seed : std::random_device{}());
	std::uniform_real_distribution<double> dst;

	std::for_each(fields.begin(), fields.end(), [&rng, &dst](double& d) {d = dst(rng); });
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

				energy += -K_s * (cos(kPi * (fields[i_r] - fields[siteIndex])) + cos(kPi * (fields[i_u] - fields[siteIndex])));
				energy += -K_t * (cos(kPi * (fields[i_t] - fields[siteIndex])));
			}
		}
	}

	return energy;
}

double QXYSquare::proposeUpdate(int index, double delta) const
{
	std::pair<std::vector<int>, std::vector<int>> connectedSites = getSiteConnectedCluster(index);

	double flip_energy = 0.0;
	for (const int& csite : connectedSites.first)
	{
		flip_energy += -K_s * (cos(kPi * (fields[index] + delta - fields[csite])) - cos(kPi * (fields[index] - fields[csite])));
	}
	for (const int& tsite : connectedSites.second)
	{
		double dTheta = fields[index] - fields[tsite];
		double dTheta_f = dTheta + delta;
		flip_energy += -K_t * (cos(kPi * dTheta_f) - cos(kPi * dTheta));
	}

	return flip_energy;
}

std::vector<std::pair<std::vector<int>, int>> QXYSquare::getSpacialVortices() const
{
	std::vector<std::pair<std::vector<int>, int>> vortices;
	for (int n = 0; n < nPlaqs; n += 3)
	{
		std::vector<int> plaq_sites = getPlaqConnectedCluster(n).first;

		double vortex = 0.0;
		size_t size = plaq_sites.size();
		for (int i = 0; i < size; i++)
		{
			int site_2 = plaq_sites[(i + 1) % size];
			int site_1 = plaq_sites[i];

			double d2 = fields[site_2];
			double d1 = fields[site_1];

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
	for (int n = 0; n < 2*nPlaqs/3; n++)
	{
		int plaq_index = n + (n - 1) / 2; // this ignores spacial plaquettes
		std::vector<int> plaq_sites = getPlaqConnectedCluster(plaq_index).first;

		double vortex = 0.0;
		size_t size = plaq_sites.size();
		for (int i = 0; i < size; i++)
		{
			int site_2 = plaq_sites[(i + 1) % size];
			int site_1 = plaq_sites[i];

			double d2 = fields[site_2];
			double d1 = fields[site_1];

			vortex += mcf::mapToCircle(d2 - d1);
		}

		if (vortex >= 1.0 - 1e-5 || vortex <= -1.0 + 1e-5)
		{
			vortices.push_back(std::pair<std::vector<int>, int>(plaq_sites, mcf::sgn(vortex)));
		}
	}

	return vortices;
}

std::vector<std::string> QXYSquare::observableNames() const
{
	return { "Energy" };
}

std::vector<double> QXYSquare::measure(double temperature) const
{
	return { getEnergy() };
}

std::vector<Channel> QXYSquare::channels() const
{
	return { { "theta", ChannelKind::Angle } };
}

void QXYSquare::fillChannel(int channel, std::vector<double>& out) const
{
	out.assign(lattice.size(), 0.0);
	for (int site = 0; site < lattice.size(); site++)
		out[site] = fields[site];
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