#include "Spiderweb.h"

#include <assert.h>

#include "MathUtil.h"
#include "SpiderwebInstantons.h"

// A fields are treated as:
// 0: A0
// 1: Axx
// 2: Axy

using mcf::kPi;

Spiderweb::Spiderweb(int linear_size, int temporal_size, double KU)
	:
	System(mcf::PeriodicLattice({ linear_size, linear_size, temporal_size }, { "x", "y", "t" }),
		linear_size * linear_size * temporal_size * 3),
	linear_size(linear_size),
	spatial_size(linear_size * linear_size),
	temporal_size(temporal_size),
	nSites(linear_size * linear_size * temporal_size),
	KU(KU)
{}

void Spiderweb::randomize(std::mt19937& rng)
{
	std::uniform_real_distribution<double> dst(-1.0, 1.0);
	std::for_each(fields.begin(), fields.end(), [&](double& d) { d = dst(rng); });
}

double Spiderweb::getEnergy() const
{
	std::vector<double> localEnergies(nSites);
	return accumulateLocalEnergies(localEnergies);
}

std::vector<double> Spiderweb::getLocalEnergies() const
{
	std::vector<double> localEnergies(nSites);
	accumulateLocalEnergies(localEnergies);
	return localEnergies;
}

double Spiderweb::getEnergy(std::vector<double>& localEnergies) const
{
	assert(localEnergies.size() == nSites);
	return accumulateLocalEnergies(localEnergies);
}

double Spiderweb::accumulateLocalEnergies(std::vector<double>& localEnergies) const
{
	// L = 1/2K cos(Q_ij A0 - D0 A_ij) - U/2 cos(Q_ij A_ij)
	// treat U as inverse temperature, such that T --> infty will make magnetic fluxes proliferate
	// Hence H = 1/2(KU) cos(Q_ij A_0 - D0 A_ij) - 1/2 cos(Q_ij A_ij)

	// Since the hamiltonian should be E^2 + B^2 and E^2 = Ex^2 + Ey^2 etc. Its probably
	// cos(Q_xx A0 - D0 A_xx) + cos(Q_xy A0 - D0 A_xy)
	double energy = 0.0;
	for (int n_site = 0; n_site < nSites; n_site++)
	{
		// (field_index, factor) pairs, such that factor * fields[field_index] is a term
		// in the cosine of the hamiltonian
		auto e_terms_xx = getElectricTerms_xx(n_site * 3);
		auto e_terms_xy = getElectricTerms_xy(n_site * 3);
		auto b_terms = getMagneticTerms(n_site * 3);

		double e_sum_xx = 0.0;
		for (auto term : e_terms_xx)
		{
			e_sum_xx += term.second * fields[term.first];
		}
		double e_sum_xy = 0.0;
		for (auto term : e_terms_xy)
		{
			e_sum_xy += term.second * fields[term.first];
		}

		double b_sum = 0.0;
		for (auto term : b_terms)
		{
			b_sum += term.second * fields[term.first];
		}

		double b_val = cos(kPi * b_sum) / 2.0;
		double local_energy = -cos(kPi * e_sum_xx) / (2.0 * KU) - cos(kPi * e_sum_xy) / (2.0 * KU) - b_val;
		energy += local_energy;

		localEnergies[n_site] = local_energy;
	}

	return energy;
}

double Spiderweb::proposeUpdate(int index, double delta) const
{
	double d_energy = 0.0;
	int type = index % 3;
	int site_index = index / 3;
	const auto site_vector = index_from_site(site_index);
	int nx = site_vector[0];
	int ny = site_vector[1];
	int nt = site_vector[2];

	int nx_m1 = lattice.wrap(0, nx - 1);
	int nx_m2 = lattice.wrap(0, nx - 2);
	int ny_m1 = lattice.wrap(1, ny - 1);
	int ny_m2 = lattice.wrap(1, ny - 2);
	int nt_m1 = lattice.wrap(2, nt - 1);

	auto accumulate_xx = [&](int anchor_field_index, double prefactor)
		{
			auto terms = getElectricTerms_xx(anchor_field_index);
			double old_sum = 0.0, new_sum = 0.0;
			for (auto term : terms)
			{
				old_sum += term.second * fields[term.first];
				new_sum += term.second * fields[term.first];
				if (term.first == index)
					new_sum += term.second * delta;
			}
			d_energy += prefactor * (-cos(kPi * new_sum) + cos(kPi * old_sum));
		};
	auto accumulate_xy = [&](int anchor_field_index, double prefactor)
		{
			auto terms = getElectricTerms_xy(anchor_field_index);
			double old_sum = 0.0, new_sum = 0.0;
			for (auto term : terms)
			{
				old_sum += term.second * fields[term.first];
				new_sum += term.second * fields[term.first];
				if (term.first == index)
					new_sum += term.second * delta;
			}
			d_energy += prefactor * (-cos(kPi * new_sum) + cos(kPi * old_sum));
		};
	auto accumulate_b = [&](int anchor_field_index, double prefactor)
		{
			auto terms = getMagneticTerms(anchor_field_index);
			double old_sum = 0.0, new_sum = 0.0;
			for (auto term : terms)
			{
				old_sum += term.second * fields[term.first];
				new_sum += term.second * fields[term.first];
				if (term.first == index)
					new_sum += term.second * delta;
			}
			d_energy += prefactor * (cos(kPi * new_sum) - cos(kPi * old_sum));
		};

	const double e_pref = 1.0 / (2.0 * KU);
	const double b_pref = -0.5;

	if (type == 0) // A_0(r): appears in E_xx at r, r-x, r-y, r-x-y; in E_xy at r-2x, r-x, r-2y, r-y
	{
		accumulate_xx(field_index_from_site(nx, ny, nt, 0), e_pref);
		accumulate_xx(field_index_from_site(nx_m1, ny, nt, 0), e_pref);
		accumulate_xx(field_index_from_site(nx, ny_m1, nt, 0), e_pref);
		accumulate_xx(field_index_from_site(nx_m1, ny_m1, nt, 0), e_pref);

		accumulate_xy(field_index_from_site(nx_m2, ny, nt, 0), e_pref);
		accumulate_xy(field_index_from_site(nx_m1, ny, nt, 0), e_pref);
		accumulate_xy(field_index_from_site(nx, ny_m2, nt, 0), e_pref);
		accumulate_xy(field_index_from_site(nx, ny_m1, nt, 0), e_pref);
	}
	else if (type == 1) // A_xx(r): appears in E_xx at r, r-t; in B at r-2x, r-x, r-2y, r-y
	{
		accumulate_xx(field_index_from_site(nx, ny, nt, 0), e_pref);
		accumulate_xx(field_index_from_site(nx, ny, nt_m1, 0), e_pref);

		accumulate_b(field_index_from_site(nx_m2, ny, nt, 0), b_pref);
		accumulate_b(field_index_from_site(nx_m1, ny, nt, 0), b_pref);
		accumulate_b(field_index_from_site(nx, ny_m2, nt, 0), b_pref);
		accumulate_b(field_index_from_site(nx, ny_m1, nt, 0), b_pref);
	}
	else // type == 2, A_xy(r): appears in E_xy at r, r-t; in B at r-x-y, r-x, r-y, r
	{
		accumulate_xy(field_index_from_site(nx, ny, nt, 0), e_pref);
		accumulate_xy(field_index_from_site(nx, ny, nt_m1, 0), e_pref);

		accumulate_b(field_index_from_site(nx_m1, ny_m1, nt, 0), b_pref);
		accumulate_b(field_index_from_site(nx_m1, ny, nt, 0), b_pref);
		accumulate_b(field_index_from_site(nx, ny_m1, nt, 0), b_pref);
		accumulate_b(field_index_from_site(nx, ny, nt, 0), b_pref);
	}

	return d_energy;
}

std::vector<std::string> Spiderweb::observableNames() const
{
	return { "Energy" };
}

std::vector<double> Spiderweb::measure(double temperature) const
{
	return { getEnergy() };
}

std::vector<Channel> Spiderweb::channels() const
{
	return { { "A_0", ChannelKind::Angle },
			 { "A_xx", ChannelKind::Angle },
			 { "A_xy", ChannelKind::Angle },
			 { "local energy", ChannelKind::Magnitude },
			 { "instantons m", ChannelKind::Integer },
			 { "Z4 steps", ChannelKind::Integer } };
}

void Spiderweb::fillChannel(int channel, std::vector<double>& out) const
{
	if (channel == 3)
	{
		out = getLocalEnergies();
		return;
	}
	// Instanton diagnostics, see SpiderwebInstantons.h. Pure measurements of the current fields.
	if (channel == 4 || channel == 5)
	{
		const std::vector<int> values = channel == 4 ? mcf::spiderwebInstantons(*this) : mcf::spiderwebZ4Steps(*this);
		out.assign(values.begin(), values.end());
		return;
	}

	out.assign(lattice.size(), 0.0);
	for (int site = 0; site < lattice.size(); site++)
		out[site] = get_field(site, channel);
}

/// <summary>
/// Returns vector of (field_index, factor) pairs, such that factor * fields[field_index] is the
/// corresponding term in the term 1/2(KU) cos(Q_xx A_0 - D0 A_xx)
/// </summary>
/// <param name="site_index"></param>
/// <returns></returns>
std::vector<std::pair<int, double>> Spiderweb::getElectricTerms_xx(int field_index) const
{
	int site_index = field_index / 3;
	const auto site_vector = index_from_site(site_index);
	int nx = site_vector[0];
	int ny = site_vector[1];
	int nt = site_vector[2];

	std::vector<std::pair<int, double>> electric_terms{};

	// -D0 A_xx = -A_xx(r + t) + A_xx(r)
	electric_terms.push_back({field_index_from_site(nx, ny, lattice.wrap(2, nt + 1), 1), -1});
	electric_terms.push_back({field_index_from_site(nx, ny, nt, 1), 1});

	// Q_xx A_00 = 4DxDy A_0 = 4A_0(r + x + y) - 4A_0(r + x) - 4A_0(r + y) + 4A_0(r)
	electric_terms.push_back({ field_index_from_site(lattice.wrap(0, nx + 1), lattice.wrap(1, ny + 1), nt, 0), 4 });
	electric_terms.push_back({ field_index_from_site(lattice.wrap(0, nx + 1), ny, nt, 0), -4 });
	electric_terms.push_back({ field_index_from_site(nx, lattice.wrap(1, ny + 1), nt, 0), -4 });
	electric_terms.push_back({ field_index_from_site(nx, ny, nt, 0), 4 });

	return electric_terms;
}

// A fields are treated as:
// 0: A0
// 1: Axx
// 2: Axy

std::vector<std::pair<int, double>> Spiderweb::getElectricTerms_xy(int field_index) const
{
	int site_index = field_index / 3;
	const auto site_vector = index_from_site(site_index);
	int nx = site_vector[0];
	int ny = site_vector[1];
	int nt = site_vector[2];

	std::vector<std::pair<int, double>> electric_terms{};

	// -D0 A_xy = -A_xy(r + t) + A_xy(r)
	electric_terms.push_back({field_index_from_site(nx, ny, lattice.wrap(2, nt + 1), 2), -1});
	electric_terms.push_back({field_index_from_site(nx, ny, nt, 2), 1});

	// Q_xy A_0 = (DxDx - DyDy)A_0 = A_0(r + 2x) - 2A_0(r + x) - A_0(r + 2y) + 2A_0(r + y)
	electric_terms.push_back({ field_index_from_site(lattice.wrap(0, nx + 2), ny, nt, 0), 1 });
	electric_terms.push_back({ field_index_from_site(lattice.wrap(0, nx + 1), ny, nt, 0), -2 });
	electric_terms.push_back({ field_index_from_site(nx, lattice.wrap(1, ny + 2), nt, 0), -1 });
	electric_terms.push_back({ field_index_from_site(nx, lattice.wrap(1, ny + 1), nt, 0), 2 });

	return electric_terms;
}

std::vector<std::pair<int, double>> Spiderweb::getMagneticTerms(int field_index) const
{
	int site_index = field_index / 3;
	const auto site_vector = index_from_site(site_index);
	int nx = site_vector[0];
	int ny = site_vector[1];
	int nt = site_vector[2];

	std::vector<std::pair<int, double>> magnetic_terms{};

	// Q_ij A_ij = (DxDx - DyDy) A_xx - 4DxDy A_xy

	// (DxDx - DyDy)A_xx = A_xx(r + 2x) - 2A_xx(r + x) - A_xx(r + 2y) + 2A_xx(r + y)
	magnetic_terms.push_back({field_index_from_site(lattice.wrap(0, nx + 2), ny, nt, 1), 1});
	magnetic_terms.push_back({field_index_from_site(lattice.wrap(0, nx + 1), ny, nt, 1), -2});
	magnetic_terms.push_back({field_index_from_site(nx, lattice.wrap(1, ny + 2), nt, 1), -1});
	magnetic_terms.push_back({field_index_from_site(nx, lattice.wrap(1, ny + 1), nt, 1), 2});

	// -4DxDy A_xy = -4A_xy(r + x + y) + 4A_xy(r + x) + 4A_xy(r + y) - 4A_xy(r)
	magnetic_terms.push_back({field_index_from_site(lattice.wrap(0, nx + 1), lattice.wrap(1, ny + 1), nt, 2), -4});
	magnetic_terms.push_back({field_index_from_site(lattice.wrap(0, nx + 1), ny, nt, 2), 4});
	magnetic_terms.push_back({field_index_from_site(nx, lattice.wrap(1, ny + 1), nt, 2), 4});
	magnetic_terms.push_back({field_index_from_site(nx, ny, nt, 2), -4});

	return magnetic_terms;
}

std::array<int, 3> Spiderweb::index_from_site(int site_index) const
{
	return std::array<int, 3>({ lattice.coord(site_index, 0),
								lattice.coord(site_index, 1),
								lattice.coord(site_index, 2) });
}

int Spiderweb::field_index_from_site(int nx, int ny, int nt, int type) const
{
	return lattice.index(nx, ny, nt) * 3 + type;
}

int Spiderweb::field_index_from_site(int site_index, int type) const
{
	return site_index * 3 + type;
}

double Spiderweb::get_field(int site_index, int type) const
{
	return fields[site_index * 3 + type];
}

double Spiderweb::get_field(int nx, int ny, int nt, int type) const
{
	return get_field(lattice.index(nx, ny, nt), type);
}