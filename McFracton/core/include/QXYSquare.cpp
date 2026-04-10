#include "QXYSquare.h"

#include <algorithm>

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
	System(size * size * Ntau, size * size * Ntau * 3)
{
	site_fields = std::vector<double>(nSites);
	plaq_fields = std::vector<double>(nSites);

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
		flip_energy +=  -K_s * (cos(2.0 * PI * (site_fields[index] + angle - site_fields[csite])) - cos(2.0 * PI * (site_fields[index] - site_fields[csite])));
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

std::vector<std::pair<std::vector<int>, int>> QXYSquare::getVortices() const
{
	return std::vector<std::pair<std::vector<int>, int>>();
}

void QXYSquare::LogToFile(std::ofstream& outfile) const
{
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

	int i_u = (nt * size + ny_u ) * size + nx;
	int i_r = (nt * size + ny   ) * size + nx_r;
	int i_d = (nt * size + ny_d ) * size + nx;
	int i_l = (nt * size + ny   ) * size + nx_l;

	int it_u = (nt_u * size + ny) * size + nx;
	int it_d = (nt_d * size + ny) * size + nx;

	return std::pair<std::vector<int>, std::vector<int>>({ i_u, i_r, i_d, i_l }, { it_u, it_d });
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