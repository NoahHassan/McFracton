#include "Spiderweb.h"

// A fields are treated as:
// 0: A0
// 1: Axx
// 2: Axy

Spiderweb::Spiderweb(int linear_size, int temporal_size, double KU)
	:
	linear_size(linear_size),
	spatial_size(linear_size * linear_size),
	temporal_size(temporal_size),
	KU(KU),
	nSites(linear_size * linear_size * temporal_size),
	nPlaqs(0),
	System(linear_size * linear_size * temporal_size * 3, 0)
{
	site_fields = std::vector<double>(n_site_variables);
	plaq_fields = std::vector<double>(n_plaq_variables);

	std::random_device rd;
	rng = std::mt19937(rd());
	overrelax_dst = std::uniform_real_distribution<double>(-1.0, 1.0);

	std::for_each(site_fields.begin(), site_fields.end(), [&](double& d) { d = overrelax_dst(rng); });
}

double Spiderweb::getEnergy() const
{
	// L = 1/2K cos(Q_ij A0 - D0 A_ij) - U/2 cos(Q_ij A_ij)
	// treat U as inverse temperature, such that T --> infty will make magnetic fluxes proliferate
	// Hence H = 1/2(KU) cos(Q_ij A_0 - D0 A_ij) - 1/2 cos(Q_ij A_ij)

	// Since the hamiltonian should be E^2 + B^2 and E^2 = Ex^2 + Ey^2 etc. Its probably
	// cos(Q_xx A0 - D0 A_xx) + cos(Q_xy A0 - D0 A_xy)
	double energy = 0.0;
	for (int n_site = 0; n_site < nSites; n_site++)
	{
		int anchor = n_site * 3;	// any field index at this site works; the term generators only use field_index/3
		auto e_terms_xx = getElectricTerms_xx(anchor);	// should return vector of (field_index, factor) pairs such that factor * site_fields[field_index]
		auto e_terms_xy = getElectricTerms_xy(anchor);	// is a term in the cosine of the hamiltonian
		auto b_terms = getMagneticTerms(anchor);

		double e_sum_xx = 0.0;
		for (auto term : e_terms_xx)
		{
			e_sum_xx += term.second * site_fields[term.first];
		}
		double e_sum_xy = 0.0;
		for (auto term : e_terms_xy)
		{
			e_sum_xy += term.second * site_fields[term.first];
		}

		double b_sum = 0.0;
		for (auto term : b_terms)
		{
			b_sum += term.second * site_fields[term.first];
		}

		energy += cos(e_sum_xx) / (2.0 * KU) + cos(e_sum_xy) / (2.0 * KU) - cos(b_sum) / 2.0;
	}

	return energy;
}

double Spiderweb::proposeSiteFlip(int index, double angle) const
{
	// IMPORTANT: because Q_ij involves second derivatives (D_x^2, D_y^2, 4 D_xD_y),
	// a given field does NOT only appear in the electric/magnetic term centered at
	// its own site r. It also appears in the terms centered at several NEIGHBORING
	// sites (the same way a term like (D_x^2 f)(s) = f(s+2x)-2f(s+x)+f(s) means that
	// changing f(r) affects the term evaluated at s=r, s=r-x AND s=r-2x).
	// Below we enumerate every (site, term-type) pair whose cosine argument actually
	// contains site_fields[index], and accumulate cos(new)-cos(old) for each.

	double d_energy = 0.0;
	int type = index % 3;
	int site_index = index / 3;
	const auto site_vector = index_from_site(site_index);
	int nx = site_vector[0];
	int ny = site_vector[1];
	int nt = site_vector[2];

	int nx_m1 = (nx - 1 + linear_size) % linear_size;
	int nx_m2 = (nx - 2 + linear_size) % linear_size;
	int ny_m1 = (ny - 1 + linear_size) % linear_size;
	int ny_m2 = (ny - 2 + linear_size) % linear_size;
	int nt_m1 = (nt - 1 + temporal_size) % temporal_size;

	auto accumulate_xx = [&](int anchor_field_index, double prefactor)
		{
			auto terms = getElectricTerms_xx(anchor_field_index);
			double old_sum = 0.0, new_sum = 0.0;
			for (auto term : terms)
			{
				old_sum += term.second * site_fields[term.first];
				new_sum += term.second * site_fields[term.first];
				if (term.first == index)
					new_sum += term.second * angle;
			}
			d_energy += prefactor * (cos(new_sum) - cos(old_sum));
		};
	auto accumulate_xy = [&](int anchor_field_index, double prefactor)
		{
			auto terms = getElectricTerms_xy(anchor_field_index);
			double old_sum = 0.0, new_sum = 0.0;
			for (auto term : terms)
			{
				old_sum += term.second * site_fields[term.first];
				new_sum += term.second * site_fields[term.first];
				if (term.first == index)
					new_sum += term.second * angle;
			}
			d_energy += prefactor * (cos(new_sum) - cos(old_sum));
		};
	auto accumulate_b = [&](int anchor_field_index, double prefactor)
		{
			auto terms = getMagneticTerms(anchor_field_index);
			double old_sum = 0.0, new_sum = 0.0;
			for (auto term : terms)
			{
				old_sum += term.second * site_fields[term.first];
				new_sum += term.second * site_fields[term.first];
				if (term.first == index)
					new_sum += term.second * angle;
			}
			d_energy += prefactor * (cos(new_sum) - cos(old_sum));
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

double Spiderweb::proposePlaqFlip(int index, double angle) const
{
	return 0.0;
}

void Spiderweb::UpdateSite(int index, double angle)
{
	site_fields[index] += angle;
}

void Spiderweb::UpdatePlaq(int index, double angle)
{
	plaq_fields[index] += angle;
}

void Spiderweb::OverrelaxSite(int index)
{
	// Do a gauge transformation A_ij --> Q_ij f with weird Q_ij
}

System::Observables Spiderweb::Measure(double T) const
{
	System::Observables observables;
	observables.energy = getEnergy();
	observables.helicity_modulus = 0.0;
	observables.n_defects_a = 0;
	observables.n_defects_b = 0;
	observables.polyakov_loop = 0.0;

	return observables;
}

/// <summary>
/// Returns vector of (field_index, factor) pairs, such that factor * site_fields[field_index] is the
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
	electric_terms.push_back({field_index_from_site(nx, ny, (nt + 1) % temporal_size, 1), -1});
	electric_terms.push_back({field_index_from_site(nx, ny, nt, 1), 1});

	// Q_xx A_00 = 4DxDy A_0 = 4A_0(r + x + y) - 4A_0(r + x) - 4A_0(r + y) + 4A_0(r)
	electric_terms.push_back({ field_index_from_site((nx + 1) % linear_size, (ny + 1) % linear_size, nt, 0), 4 });
	electric_terms.push_back({ field_index_from_site((nx + 1) % linear_size, ny, nt, 0), -4 });
	electric_terms.push_back({ field_index_from_site(nx, (ny + 1) % linear_size, nt, 0), -4 });
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
	electric_terms.push_back({field_index_from_site(nx, ny, (nt + 1) % temporal_size, 2), -1});
	electric_terms.push_back({field_index_from_site(nx, ny, nt, 2), 1});

	// Q_xy A_0 = (DxDx - DyDy)A_0 = A_0(r + 2x) - 2A_0(r + x) - A_0(r + 2y) + 2A_0(r + y)
	electric_terms.push_back({ field_index_from_site((nx + 2) % linear_size, ny, nt, 0), 1 });
	electric_terms.push_back({ field_index_from_site((nx + 1) % linear_size, ny, nt, 0), -2 });
	electric_terms.push_back({ field_index_from_site(nx, (ny + 2) % linear_size, nt, 0), -1 });
	electric_terms.push_back({ field_index_from_site(nx, (ny + 1) % linear_size, nt, 0), 2 });

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
	magnetic_terms.push_back({field_index_from_site((nx + 2) % linear_size, ny, nt, 1), 1});
	magnetic_terms.push_back({field_index_from_site((nx + 1) % linear_size, ny, nt, 1), -2});
	magnetic_terms.push_back({field_index_from_site(nx, (ny + 2) % linear_size, nt, 1), -1});
	magnetic_terms.push_back({field_index_from_site(nx, (ny + 1) % linear_size, nt, 1), 2});

	// -4DxDy A_xy = -4A_xy(r + x + y) + 4A_xy(r + x) + 4A_xy(r + y) - 4A_xy(r)
	magnetic_terms.push_back({field_index_from_site((nx + 1) % linear_size, (ny + 1) % linear_size, nt, 2), -4});
	magnetic_terms.push_back({field_index_from_site((nx + 1) % linear_size, ny, nt, 2), 4});
	magnetic_terms.push_back({field_index_from_site(nx, (ny + 1) % linear_size, nt, 2), 4});
	magnetic_terms.push_back({field_index_from_site(nx, ny, nt, 2), -4});

	return magnetic_terms;
}

const int Spiderweb::to_site_index(int nx, int ny, int nt) const
{
	return (nt * linear_size + ny) * linear_size + nx;
}

const std::array<int, 3> Spiderweb::index_from_site(int site_index) const
{
	int nt = site_index / spatial_size;
	int ny = (site_index % spatial_size) / linear_size;
	int nx = site_index % linear_size;

	return std::array<int, 3>({ nx, ny, nt });
}

const int Spiderweb::field_index_from_site(int nx, int ny, int nt, int type) const
{
	return to_site_index(nx, ny, nt) * 3 + type;
}

const int Spiderweb::field_index_from_site(int site_index, int type) const
{
	return site_index * 3 + type;
}

double Spiderweb::get_field(int site_index, int type) const
{
	return site_fields[site_index + type * 3];
}

double Spiderweb::get_field(int nx, int ny, int nt, int type) const
{
	return get_field(to_site_index(nx, ny, nt), type);
}

double Spiderweb::mapToCircle(const double& d) const
{
	if (d >= 0.5)
		return d - int(d + 0.5);
	else if (d <= 0.5)
		return d - int(d - 0.5);
	else
		return d;
}