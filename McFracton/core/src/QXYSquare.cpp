#include "QXYSquare.h"

#include <algorithm>
#include <numeric>
#include <assert.h>

#ifndef PI
#define PI 3.1415926535897932384
#endif

QXYSquare::QXYSquare(int size, int Ntau)
	:
	QXYSquare(size, Ntau, 1.0, 1.0)
{}

QXYSquare::QXYSquare(int size, int Ntau, float K_s, float K_t)
	:
	size(size),
	Ntau(Ntau),
	K_s(K_s),
	K_t(K_t),
	ss_size(size * size),
	st_size(size * Ntau),
	System(size * size * Ntau, size * size * Ntau * 3)
{
	site_fields = std::vector<double>(nSites);
	plaq_fields = std::vector<double>(nSites*3);

	std::random_device rd;
	std::mt19937 rng(rd());
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
				int siteIndex = (nt * size + ny) * size + nx;

				// No double counting
				int n_r = (nx + 1) % size;
				int n_u = (ny + 1) % size;
				int n_t = (nt + 1) % Ntau;

				int i_r = (nt * size + ny) * size + n_r;
				int i_u = (nt * size + n_u) * size + nx;
				int i_t = (n_t * size + ny) * size + nx;

				energy += -K_s * (cos(2.0 * PI * (site_fields[i_r] - site_fields[siteIndex])) + cos(2.0 * PI * (site_fields[i_u] - site_fields[siteIndex])));
				energy += -K_t * (cos(2.0 * PI * (site_fields[i_t] - site_fields[siteIndex])));
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
		flip_energy += -K_s * (cos(2.0 * PI * (site_fields[index] + angle - site_fields[csite])) - cos(2.0 * PI * (site_fields[index] - site_fields[csite])));
	}
	for (const int& tsite : connectedSites.second)
	{
		double dTheta = site_fields[index] - site_fields[tsite];
		double dTheta_f = dTheta + angle;
		flip_energy += -K_t * (cos(2.0 * PI * dTheta_f) - cos(2.0 * PI * dTheta));
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
	for (int n = 0; n < nPlaqs; n += 3)
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

			vortex += mapToCircle(d2 - d1);
		}

		if (vortex >= 1.0 - 1e-5 || vortex <= -1.0 + 1e-5)
		{
			vortices.push_back(std::pair<std::vector<int>, int>(plaq_sites, sgn(vortex)));
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

			double d2 = site_fields[site_2];
			double d1 = site_fields[site_1];

			vortex += mapToCircle(d2 - d1);
		}

		if (vortex >= 1.0 - 1e-5 || vortex <= -1.0 + 1e-5)
		{
			vortices.push_back(std::pair<std::vector<int>, int>(plaq_sites, sgn(vortex)));
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

const std::pair<std::vector<int>, std::vector<int>> QXYSquare::getSiteConnectedCluster(int siteIndex) const
{
	int nt = siteIndex / (size * size);
	int ny = (siteIndex - nt * size * size) / size;
	int nx = siteIndex - nt * size * size - ny * size;

	int nx_r = (nx + 1) % size;
	int nx_l = (nx - 1 + size) % size;
	int ny_u = (ny + 1) % size;
	int ny_d = (ny - 1 + size) % size;

	int nt_u = (nt + 1) % Ntau;
	int nt_d = (nt - 1 + Ntau) % Ntau;

	int i_u = (nt * size + ny_u) * size + nx;
	int i_r = (nt * size + ny) * size + nx_r;
	int i_d = (nt * size + ny_d) * size + nx;
	int i_l = (nt * size + ny) * size + nx_l;

	int it_u = (nt_u * size + ny) * size + nx;
	int it_d = (nt_d * size + ny) * size + nx;

	return std::pair<std::vector<int>, std::vector<int>>({ i_u, i_r, i_d, i_l }, { it_u, it_d });
}

const std::pair<std::vector<int>, std::vector<int>> QXYSquare::getPlaqConnectedCluster(int plaqIndex) const
{
	int plaqType = plaqIndex % 3;
	int plaqSite = plaqIndex / 3;

	int n_x = plaqSite % size;
	int n_y = (plaqSite % ss_size) / size;
	int n_t = plaqSite / ss_size;

	int i_a = plaqSite;
	int i_b = -1;
	int i_c = -1;
	int i_d = -1;

	auto to_index = [&](int mx, int my, int mz) {
		assert(mz < Ntau);
		assert(my < size);
		assert(mx < size);
		return mz * ss_size + my * size + mx;
		};

	switch (plaqType)
	{
	case 0:
		i_b = to_index(n_x, (n_y + 1) % size, n_t);
		i_c = to_index((n_x + 1) % size, (n_y + 1) % size, n_t);
		i_d = to_index((n_x + 1) % size, n_y, n_t);
		break;
	case 1:
		i_b = to_index(n_x, n_y, (n_t + 1) % Ntau);
		i_c = to_index((n_x + 1) % size, n_y, (n_t + 1) % Ntau);
		i_d = to_index((n_x + 1) % size, n_y, n_t);
		break;
	case 2:
		i_b = to_index(n_x, n_y, (n_t + 1) % Ntau);
		i_c = to_index(n_x, (n_y + 1) % size, (n_t + 1) % Ntau);
		i_d = to_index(n_x, (n_y + 1) % size, n_t);
		break;
	}

	return std::pair<std::vector<int>, std::vector<int>>({ i_a, i_b, i_c, i_d }, {});
}

double QXYSquare::mapToCircle(const double& d) const
{
	if (d >= 0.5)
		return d - int(d + 0.5);
	else if (d <= 0.5)
		return d - int(d - 0.5);
	else
		return d;
}