#include "XYSquare.h"

#include <cmath>
#include <assert.h>

#ifndef PI
#define PI 3.1415926535897932384
#endif

XYSquare::XYSquare(int size)
	:
	System(size * size, (size-1) * (size-1)),
	length(size)
{
	site_fields = std::vector<double>(nSites);
	plaq_fields = std::vector<double>(nSites);
}

double XYSquare::getEnergy() const
{
	double energy = 0.0;
	for (int nx = 0; nx < length; nx++)
	{
		for (int ny = 0; ny < length; ny++)
		{
			int siteIndex = ny * length + nx;

			int nx_r = (nx + 1) % length;
			int ny_u = (ny + 1) % length;

			int i_r = ny * length + nx_r;
			int i_u = ny_u * length + nx;

			energy += cos(2.0 * PI * (site_fields[i_r] - site_fields[siteIndex])) + cos(2.0 * PI * (site_fields[i_u] - site_fields[siteIndex]));
		}
	}

	return -energy;
}

double XYSquare::proposeSiteFlip(int index, double angle) const
{
	std::vector<int> connectedSites = getSiteConnectedCluster(index).first;

	double flip_energy = 0.0;
	for (const int& csite : connectedSites) {
		flip_energy += cos(2.0 * PI * (site_fields[index] + angle - site_fields[csite])) - cos(2.0 * PI * (site_fields[index] - site_fields[csite]));
	}

	// flip plaquettes too

	return -flip_energy;
}

double XYSquare::proposePlaqFlip(int index, double angle) const
{
	return 0.0;
}

void XYSquare::UpdateSite(int index, double angle)
{
	site_fields[index] += angle;
}

void XYSquare::UpdatePlaq(int index, double angle)
{
	plaq_fields[index] += angle;
}

std::vector<std::pair<std::vector<int>, int>> XYSquare::getVortices() const
{
	std::vector<std::pair<std::vector<int>, int>> vortices;
	for (int n = 0; n < nPlaqs; n++)
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

const std::pair<std::vector<int>, std::vector<int>> XYSquare::getSiteConnectedCluster(int siteIndex) const
{
	int ny = siteIndex / length;
	int nx = siteIndex - ny * length;

	int nx_r = (nx + 1) % length;
	int nx_l = (nx - 1 + length) % length;
	int ny_u = (ny + 1) % length;
	int ny_d = (ny - 1 + length) % length;

	int i_u = ny_u * length + nx;
	int i_r = ny * length + nx_r;
	int i_d = ny_d * length + nx;
	int i_l = ny * length + nx_l;
	return std::pair<std::vector<int>, std::vector<int>>({ i_u, i_r, i_d, i_l }, {});
}

const std::pair<std::vector<int>, std::vector<int>> XYSquare::getPlaqConnectedCluster(int plaqIndex) const
{
	int ny = plaqIndex / (length - 1);
	int nx = plaqIndex - ny * (length - 1);

	int site_bl = length * ny + nx;
	int site_br = length * ny + nx + 1;
	int site_tl = length * (ny + 1) + nx;
	int site_tr = length * (ny + 1) + nx + 1;
	return std::pair<std::vector<int>, std::vector<int>>({ site_bl, site_br, site_tr, site_tl }, {});
}

double XYSquare::mapToCircle(const double& d) const
{
	if (d >= 0.5)
		return d - int(d + 0.5);
	else if (d <= 0.5)
		return d - int(d - 0.5);
	else
		return d;
}