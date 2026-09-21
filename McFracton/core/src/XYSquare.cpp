#include "XYSquare.h"

#include <cmath>
#include <assert.h>

#include "MathUtil.h"

using mcf::kPi;

XYSquare::XYSquare(int size)
	:
	size(size),
	lattice({ size, size }, { "x", "y" }),
	System(size * size, size * size)
{
	site_fields = std::vector<double>(n_site_variables);
	plaq_fields = std::vector<double>(n_site_variables);
}

double XYSquare::getEnergy() const
{
	double energy = 0.0;
	for (int nx = 0; nx < size; nx++)
	{
		for (int ny = 0; ny < size; ny++)
		{
			int siteIndex = lattice.index(nx, ny);

			// No double counting
			int i_r = lattice.neighbor(siteIndex, 0, +1);
			int i_u = lattice.neighbor(siteIndex, 1, +1);

			energy += cos(kPi * (site_fields[i_r] - site_fields[siteIndex])) + cos(kPi * (site_fields[i_u] - site_fields[siteIndex]));
		}
	}

	return -energy;
}

double XYSquare::getSinSqrX() const
{
	double result = 0.0;
	for (int nx = 0; nx < size; nx++)
	{
		for (int ny = 0; ny < size; ny++)
		{
			int siteIndex = lattice.index(nx, ny);

			// No double counting
			int i_r = lattice.neighbor(siteIndex, 0, +1);

			result += sin(kPi * (site_fields[i_r] - site_fields[siteIndex]));
		}
	}

	return result * result;
}

double XYSquare::proposeSiteFlip(int index, double angle) const
{
	std::vector<int> connectedSites = getSiteConnectedCluster(index).first;

	double flip_energy = 0.0;
	for (const int& csite : connectedSites) {
		flip_energy += cos(kPi * (site_fields[index] + angle - site_fields[csite])) - cos(kPi * (site_fields[index] - site_fields[csite]));
	}

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
	for (int n = 0; n < n_plaq_variables; n++)
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

		if (vortex >= 2.0 - 1e-5 || vortex <= -2.0 + 1e-5)
		{
			vortices.push_back(std::pair<std::vector<int>, int>(plaq_sites, mcf::sgn(vortex)));
		}
	}

	return vortices;
}

System::Observables XYSquare::Measure(double T) const
{
	System::Observables observables;
	observables.energy = getEnergy();
	observables.flux_cos = 0.0;
	observables.helicity_modulus = 
		-observables.energy / (2.0 * (double)n_site_variables) - 
		getSinSqrX() / (T * (double)n_site_variables);
	observables.polyakov_loop = 0.0;

	int n_a = 0;
	int n_b = 0;
	const auto monopoles = getVortices();
	for (int n = 0; n < monopoles.size(); n++)
	{
		if (monopoles[n].second < 0)
			n_b++;
		else if (monopoles[n].second > 0)
			n_a++;
	}

	observables.n_defects_a = n_a;
	observables.n_defects_b = n_b;

	return observables;
}

std::pair<std::vector<int>, std::vector<int>> XYSquare::getSiteConnectedCluster(int siteIndex) const
{
	int i_u = lattice.neighbor(siteIndex, 1, +1);
	int i_r = lattice.neighbor(siteIndex, 0, +1);
	int i_d = lattice.neighbor(siteIndex, 1, -1);
	int i_l = lattice.neighbor(siteIndex, 0, -1);
	return std::pair<std::vector<int>, std::vector<int>>({ i_u, i_r, i_d, i_l }, {});
}

// The plaquette anchored at a site, as its four corner sites counter-clockwise.
std::pair<std::vector<int>, std::vector<int>> XYSquare::getPlaqConnectedCluster(int plaqIndex) const
{
	int site_bl = plaqIndex;
	int site_br = lattice.neighbor(site_bl, 0, +1);
	int site_tl = lattice.neighbor(site_bl, 1, +1);
	int site_tr = lattice.neighbor(site_br, 1, +1);
	return std::pair<std::vector<int>, std::vector<int>>({ site_bl, site_br, site_tr, site_tl }, {});
}