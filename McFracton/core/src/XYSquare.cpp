#include "XYSquare.h"

#include <cmath>
#include <assert.h>

#include "MathUtil.h"

using mcf::kPi;

XYSquare::XYSquare(int size)
	:
	System(mcf::PeriodicLattice({ size, size }, { "x", "y" }), size * size),
	size(size)
{}

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

			energy += cos(kPi * (fields[i_r] - fields[siteIndex])) + cos(kPi * (fields[i_u] - fields[siteIndex]));
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

			result += sin(kPi * (fields[i_r] - fields[siteIndex]));
		}
	}

	return result * result;
}

double XYSquare::proposeUpdate(int index, double delta) const
{
	std::vector<int> connectedSites = getSiteConnectedCluster(index).first;

	double flip_energy = 0.0;
	for (const int& csite : connectedSites) {
		flip_energy += cos(kPi * (fields[index] + delta - fields[csite])) - cos(kPi * (fields[index] - fields[csite]));
	}

	return -flip_energy;
}

std::vector<std::pair<std::vector<int>, int>> XYSquare::getVortices() const
{
	std::vector<std::pair<std::vector<int>, int>> vortices;
	for (int n = 0; n < numVariables(); n++)
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

		if (vortex >= 2.0 - 1e-5 || vortex <= -2.0 + 1e-5)
		{
			vortices.push_back(std::pair<std::vector<int>, int>(plaq_sites, mcf::sgn(vortex)));
		}
	}

	return vortices;
}

std::vector<std::string> XYSquare::observableNames() const
{
	return { "Energy", "Helicity Modulus", "defects_a", "defects_b" };
}

std::vector<double> XYSquare::measure(double temperature) const
{
	const double energy = getEnergy();
	const double helicity_modulus =
		-energy / (2.0 * (double)numVariables()) -
		getSinSqrX() / (temperature * (double)numVariables());

	int n_a = 0;
	int n_b = 0;
	const auto vortices = getVortices();
	for (int n = 0; n < vortices.size(); n++)
	{
		if (vortices[n].second < 0)
			n_b++;
		else if (vortices[n].second > 0)
			n_a++;
	}

	return { energy, helicity_modulus, (double)n_a, (double)n_b };
}

std::vector<Channel> XYSquare::channels() const
{
	return { { "theta", ChannelKind::Angle }, { "vortices", ChannelKind::Integer } };
}

void XYSquare::fillChannel(int channel, std::vector<double>& out) const
{
	out.assign(lattice.size(), 0.0);

	if (channel == 0)
	{
		for (int site = 0; site < lattice.size(); site++)
			out[site] = fields[site];
		return;
	}

	// One charge per plaquette, written to the site the plaquette is anchored at.
	for (const auto& vortex : getVortices())
		out[vortex.first[0]] = (double)vortex.second;
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