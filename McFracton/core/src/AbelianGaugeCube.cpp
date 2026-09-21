#include "AbelianGaugeCube.h"

#include <assert.h>
#include <iostream>
#include <stdexcept>

#include "MathUtil.h"

using mcf::kPi;

AbelianGaugeCube::AbelianGaugeCube(int linear_size, int temporal_size, unsigned int seed)
	:
	linear_size(linear_size),
	temporal_size(temporal_size),
	lattice({ linear_size, linear_size, linear_size, temporal_size }, { "x", "y", "z", "t" }),
	nSites(linear_size * linear_size * linear_size * temporal_size),
	nPlaqs(linear_size * linear_size * linear_size * temporal_size * 6), // xy xz xt yz yt zt
	System(linear_size * linear_size * linear_size * temporal_size * 4, 0) // A_x, A_y, A_z, A_t on each site
{
	site_fields = std::vector<double>(n_site_variables);
	plaq_fields = std::vector<double>(n_site_variables * 3);

	// seed == 0 means "pick a fresh, unpredictable seed"; any other value is reproducible.
	rng = std::mt19937(seed != 0 ? seed : std::random_device{}());
	overrelax_dst = std::uniform_real_distribution<double>(-1.0, 1.0);

	std::for_each(site_fields.begin(), site_fields.end(), [&](double& d) { d = overrelax_dst(rng); });
}

double AbelianGaugeCube::getEnergy() const
{
	double energy = 0.0;
	for (int n_plaq = 0; n_plaq < nPlaqs; n_plaq++)
	{
		auto connected_fields = getPlaqConnectedFields(n_plaq);
		double plaquette_sum = sum_plaquette(connected_fields);
		energy += cos(kPi * plaquette_sum);
	}

	return -energy;
}

double AbelianGaugeCube::proposeSiteFlip(int index, double angle) const
{
	int site_index = index / 4;
	int type = index % 4;

	int nt = lattice.coord(site_index, 3);
	int nz = lattice.coord(site_index, 2);
	int ny = lattice.coord(site_index, 1);
	int nx = lattice.coord(site_index, 0);

	double current = site_fields[index];

	switch (type)
	{
	case 0:
	{
		return getLocalEnergy_x(nx, ny, nz, nt, current + angle) - getLocalEnergy_x(nx, ny, nz, nt, current);
		break;
	}
	case 1:
	{
		return getLocalEnergy_y(nx, ny, nz, nt, current + angle) - getLocalEnergy_y(nx, ny, nz, nt, current);
		break;
	}
	case 2:
	{
		return getLocalEnergy_z(nx, ny, nz, nt, current + angle) - getLocalEnergy_z(nx, ny, nz, nt, current);
		break;
	}
	case 3:
		return getLocalEnergy_t(nx, ny, nz, nt, current + angle) - getLocalEnergy_t(nx, ny, nz, nt, current);
		break;
	default:
		throw std::logic_error("AbelianGaugeCube: link direction must be 0, 1, 2 or 3");
	}
}

double AbelianGaugeCube::proposePlaqFlip(int index, double angle) const
{
	return plaq_fields[index];
}

void AbelianGaugeCube::UpdateSite(int index, double angle)
{
	site_fields[index] += angle;
}

void AbelianGaugeCube::UpdatePlaq(int index, double angle)
{
	plaq_fields[index] += angle;
}

void AbelianGaugeCube::OverrelaxSite(int index)
{
	int site_index = index / 4;

	int nt = lattice.coord(site_index, 3);
	int nz = lattice.coord(site_index, 2);
	int ny = lattice.coord(site_index, 1);
	int nx = lattice.coord(site_index, 0);

	double random_shift = overrelax_dst(rng);
	// A_i(r) --> A_i(r) + f(r + i) - f(r) = A_i(r) - f(r)
	site_fields[to_site_index(nx, ny, nz, nt) * 4 + 0] += -random_shift;
	site_fields[to_site_index(nx, ny, nz, nt) * 4 + 1] += -random_shift;
	site_fields[to_site_index(nx, ny, nz, nt) * 4 + 2] += -random_shift;
	site_fields[to_site_index(nx, ny, nz, nt) * 4 + 3] += -random_shift;

	// A_x-1(r) --> A_x-1(r) + f(r) - f(r-x) = A_x-1(r) + f(r)
	site_fields[to_site_index(lattice.wrap(0, nx - 1), ny, nz, nt) * 4 + 0] += random_shift;
	// A_y-1(r) --> A_y-1(r) + f(r) - f(r-y) = A_y-1(r) + f(r)
	site_fields[to_site_index(nx, lattice.wrap(1, ny - 1), nz, nt) * 4 + 1] += random_shift;
	// A_z-1(r) --> A_z-1(r) + f(r) - f(r-z) = A_z-1(r) + f(r)
	site_fields[to_site_index(nx, ny, lattice.wrap(2, nz - 1), nt) * 4 + 2] += random_shift;
	// A_t-1(r) --> A_t-1(r) + f(r) - f(r-t) = A_t-1(r) + f(r)
	site_fields[to_site_index(nx, ny, nz, lattice.wrap(3, nt - 1)) * 4 + 3] += random_shift;
}

std::vector<double> AbelianGaugeCube::getFluxes_z() const
{
	std::vector<double> fluxes(nSites);
	const int cube_size = linear_size * linear_size * linear_size;
	const int plane_size = linear_size * linear_size;
	for (int n_space = 0; n_space < cube_size; n_space++)
	{
		int nz = n_space / plane_size;
		int ny = (n_space % plane_size) / linear_size;
		int nx = n_space % linear_size;
		for (int nt = 0; nt < temporal_size; nt++)
		{
			// Because every site corresponds to 6 plaquettes
			double flux = mcf::mapToCircle(
				sum_plaquette(
					getPlaqConnectedFields(
						to_site_index(nx, ny, nz, nt) * 6 + 0
					)
				)
			);

			fluxes[to_site_index(nx, ny, nz, nt)] = flux;
		}
	}

	return fluxes;
}

System::Observables AbelianGaugeCube::Measure(double T) const
{
	System::Observables observables;
	observables.energy = getEnergy();
	observables.helicity_modulus = 0.0;

	// Monopole counting was removed with getMonopoles(); the 3+1D defect density is not measured.
	observables.n_defects_a = 0;
	observables.n_defects_b = 0;

	int plane_size = linear_size * linear_size;
	int cube_size = plane_size * linear_size;
	std::vector<double> loops(cube_size);
	for (int n_space = 0; n_space < cube_size; n_space++)
	{
		double field_sum = 0.0;
		for (int nt = 0; nt < temporal_size; nt++)
		{
			int site_index = cube_size * nt + n_space;
			field_sum += get_field(site_index, 3);
		}
		loops[n_space] = cos(kPi * field_sum);
	}

	double polyakov_loop_mean = 0.0;
	std::for_each(loops.begin(), loops.end(), [&polyakov_loop_mean](double& d) {polyakov_loop_mean += d; });
	polyakov_loop_mean /= (double)cube_size;

	observables.polyakov_loop = polyakov_loop_mean;

	return observables;
};

double AbelianGaugeCube::getLocalEnergy_x(int nx, int ny, int nz, int nt, double angle) const
{
	// 0: xy, 1: xz, 2: xt, 3: yz, 4: yt, 5: zt
	double energy = 0.0;
	std::vector<std::vector<std::pair<int, int>>> connected_plaquettes = {
		getPlaqConnectedFields(nx, ny, nz, nt, 0), // xy
		getPlaqConnectedFields(nx, ny, nz, nt, 1), // xz
		getPlaqConnectedFields(nx, ny, nz, nt, 2), // xt
		getPlaqConnectedFields(nx, lattice.wrap(1, ny - 1), nz, nt, 0),	// xy
		getPlaqConnectedFields(nx, ny, lattice.wrap(2, nz - 1), nt, 1),	// xz
		getPlaqConnectedFields(nx, ny, nz, lattice.wrap(3, nt - 1), 2) // xt
	};

	for (int n = 0; n < connected_plaquettes.size(); n++)
	{
		auto plaquette = connected_plaquettes[n];
		int angle_index = (n / 3) * 2; // because the links appear at first or third position
		double plaquette_sum = sum_plaquette(plaquette, angle_index, angle);

		energy += cos(kPi * plaquette_sum);
	}

	return -energy;
}

double AbelianGaugeCube::getLocalEnergy_y(int nx, int ny, int nz, int nt, double angle) const
{
	// 0: xy, 1: xz, 2: xt, 3: yz, 4: yt, 5: zt
	double energy = 0.0;
	std::vector<std::vector<std::pair<int, int>>> connected_plaquettes = {
		getPlaqConnectedFields(nx, ny, nz, nt, 0),										// xy
		getPlaqConnectedFields(lattice.wrap(0, nx - 1), ny, nz, nt, 0),	// xy
		getPlaqConnectedFields(nx, ny, nz, nt, 3),										// yz
		getPlaqConnectedFields(nx, ny, lattice.wrap(2, nz - 1), nt, 3),	// yz
		getPlaqConnectedFields(nx, ny, nz, nt, 4),										// yt
		getPlaqConnectedFields(nx, ny, nz, lattice.wrap(3, nt - 1), 4)	// yt
	};

	int angle_indices[] = { 3, 1, 0, 2, 0, 2 };
	for (int n = 0; n < connected_plaquettes.size(); n++)
	{
		auto plaquette = connected_plaquettes[n];
		int angle_index = angle_indices[n];
		double plaquette_sum = sum_plaquette(plaquette, angle_index, angle);

		energy += cos(kPi * plaquette_sum);
	}

	return -energy;
}

double AbelianGaugeCube::getLocalEnergy_z(int nx, int ny, int nz, int nt, double angle) const
{
	// 0: xy, 1: xz, 2: xt, 3: yz, 4: yt, 5: zt
	double energy = 0.0;
	std::vector<std::vector<std::pair<int, int>>> connected_plaquettes = {
		getPlaqConnectedFields(nx, ny, nz, nt, 1),										// xz
		getPlaqConnectedFields(lattice.wrap(0, nx - 1), ny, nz, nt, 1),	// xz
		getPlaqConnectedFields(nx, ny, nz, nt, 3),										// yz
		getPlaqConnectedFields(nx, lattice.wrap(1, ny - 1), nz, nt, 3),	// yz
		getPlaqConnectedFields(nx, ny, nz, nt, 5),										// zt
		getPlaqConnectedFields(nx, ny, nz, lattice.wrap(3, nt - 1), 5)	// zt
	};

	int angle_indices[] = { 3, 1, 3, 1, 0, 2 };
	for (int n = 0; n < connected_plaquettes.size(); n++)
	{
		auto plaquette = connected_plaquettes[n];
		int angle_index = angle_indices[n];
		double plaquette_sum = sum_plaquette(plaquette, angle_index, angle);

		energy += cos(kPi * plaquette_sum);
	}

	return -energy;
}

double AbelianGaugeCube::getLocalEnergy_t(int nx, int ny, int nz, int nt, double angle) const
{
	// 0: xy, 1: xz, 2: xt, 3: yz, 4: yt, 5: zt
	double energy = 0.0;
	std::vector<std::vector<std::pair<int, int>>> connected_plaquettes = {
		getPlaqConnectedFields(nx, ny, nz, nt, 2),										// xt
		getPlaqConnectedFields(lattice.wrap(0, nx - 1), ny, nz, nt, 2),	// xt
		getPlaqConnectedFields(nx, ny, nz, nt, 4),										// yt
		getPlaqConnectedFields(nx, lattice.wrap(1, ny - 1), nz, nt, 4),	// yt
		getPlaqConnectedFields(nx, ny, nz, nt, 5),										// zt
		getPlaqConnectedFields(nx, ny, lattice.wrap(2, nz - 1), nt, 5)		// zt
	};

	int angle_indices[] = { 3, 1, 3, 1, 3, 1 };
	for (int n = 0; n < connected_plaquettes.size(); n++)
	{
		auto plaquette = connected_plaquettes[n];
		int angle_index = angle_indices[n];
		double plaquette_sum = sum_plaquette(plaquette, angle_index, angle);

		energy += cos(kPi * plaquette_sum);
	}

	return -energy;
}

/// <summary>
/// Given a plaquette, returns the corresponding fields as pairs of site and direction.
/// </summary>
/// <param name="nx"></param>
/// <param name="ny"></param>
/// <param name="nt"></param>
/// <param name="type">0:xy, 1:tx, 2:yt</param>
/// <returns></returns>
std::vector<std::pair<int, int>> AbelianGaugeCube::getPlaqConnectedFields(int nx, int ny, int nz, int nt, int type) const
{
	// 0:xy, 1:xz, 2:xt, 3:yz, 4:yt, 5:zt
	switch (type)
	{
	case 0: // xy-plaquette: A_x(r)+A_y(r+x)-A_x(r+y)-A_y(r)
		return {
			{to_site_index(nx, ny, nz, nt), 0},							// A_x(r)
			{to_site_index(lattice.wrap(0, nx + 1), ny, nz, nt), 1},		// A_y(r+x)
			{to_site_index(nx, lattice.wrap(1, ny + 1), nz, nt), 0},		// A_x(r+y)
			{to_site_index(nx, ny, nz, nt), 1}							// A_y(r)
		};
		break;
	case 1: // xz-plaquette: A_x(r)+A_z(r+x)-A_x(r+z)-A_z(r)
		return {
			{to_site_index(nx, ny, nz, nt), 0},							// A_x(r)
			{to_site_index(lattice.wrap(0, nx + 1), ny, nz, nt), 2},		// A_z(r+x)
			{to_site_index(nx, ny, lattice.wrap(2, nz + 1), nt), 0},		// A_x(r+z)
			{to_site_index(nx, ny, nz, nt), 2}							// A_z(r)
		};
		break;
	case 2: // xt-plaquette: A_x(r)+A_t(r+x)-A_x(r+t)-A_t(r)
		return {
			{to_site_index(nx, ny, nz, nt), 0},							// A_x(r)
			{to_site_index(lattice.wrap(0, nx + 1), ny, nz, nt), 3},		// A_t(r+x)
			{to_site_index(nx, ny, nz, lattice.wrap(3, nt + 1)), 0},	// A_x(r+t)
			{to_site_index(nx, ny, nz, nt), 3}							// A_t(r)
		};
		break;
	case 3: // yz-plaquette: A_y(r)+A_z(r+y)-A_y(r+z)-A_z(r)
		return {
			{to_site_index(nx, ny, nz, nt), 1},							// A_y(r)
			{to_site_index(nx, lattice.wrap(1, ny + 1), nz, nt), 2},		// A_z(r+y)
			{to_site_index(nx, ny, lattice.wrap(2, nz + 1), nt), 1},		// A_x(r+z)
			{to_site_index(nx, ny, nz, nt), 2}							// A_z(r)
		};
		break;
	case 4: // yt-plaquette: A_y(r)+A_t(r+y)-A_y(r+t)-A_t(r)
		return {
			{to_site_index(nx, ny, nz, nt), 1},							// A_y(r)
			{to_site_index(nx, lattice.wrap(1, ny + 1), nz, nt), 3},		// A_t(r+y)
			{to_site_index(nx, ny, nz, lattice.wrap(3, nt + 1)), 1},	// A_y(r+t)
			{to_site_index(nx, ny, nz, nt), 3}							// A_t(r)
		};
		break;
	case 5: // zt-plaquette: A_z(r)+A_t(r+z)-A_z(r+t)-A_t(r)
		return {
			{to_site_index(nx, ny, nz, nt), 2},							// A_z(r)
			{to_site_index(nx, ny, lattice.wrap(2, nz + 1), nt), 3},		// A_t(r+z)
			{to_site_index(nx, ny, nz, lattice.wrap(3, nt + 1)), 2},	// A_z(r+t)
			{to_site_index(nx, ny, nz, nt), 3}							// A_t(r)
		};
		break;
	default:
		throw std::logic_error("AbelianGaugeCube: plaquette type must be 0..5");
	}
}

/// <summary>
/// Given a plaquette, returns the corresponding fields as pairs of site and direction.
/// </summary>
/// <param name="plaqIndex"></param>
/// <returns></returns>
std::vector<std::pair<int, int>> AbelianGaugeCube::getPlaqConnectedFields(int plaqIndex) const
{
	// Each site connects uniquely to six plaquettes
	// Hence plaqIndex / 6 is the site and plaqIndex % 6 is the type
	int site_index = plaqIndex / 6;
	int plaq_type = plaqIndex % 6;

	return getPlaqConnectedFields(lattice.coord(site_index, 0),
		lattice.coord(site_index, 1),
		lattice.coord(site_index, 2),
		lattice.coord(site_index, 3),
		plaq_type);
}

// e.g. A_x(r) + A_y(r+x) - A_x(r+y) - A_y(r)
double AbelianGaugeCube::sum_plaquette(const std::vector<std::pair<int, int>>& plaquette) const
{
	double sum = 0.0;
	for (int i = 0; i < plaquette.size(); i++)
	{
		double field_sign = double(1 - 2 * (i / 2)); // +1, +1, -1, -1
		sum += field_sign * get_field(plaquette[i].first, plaquette[i].second);
	}
	return sum;
}

double AbelianGaugeCube::sum_plaquette(const std::vector<std::pair<int, int>>& plaquette, int angle_index, double angle) const
{
	double sum = 0.0;
	for (int i = 0; i < plaquette.size(); i++)
	{
		double field_sign = double(1 - 2 * (i / 2)); // +1, +1, -1, -1
		if (i == angle_index)
		{
			sum += field_sign * angle;
		}
		else
		{
			sum += field_sign * get_field(plaquette[i].first, plaquette[i].second);
		}
	}
	return sum;
}

int AbelianGaugeCube::to_site_index(int nx, int ny, int nz, int nt) const
{
	return lattice.index(nx, ny, nz, nt);
}

double AbelianGaugeCube::get_field(int site_index, int direction) const
{
	assert(0 <= direction && direction <= 4);
	assert(site_index <= nSites);

	return site_fields[site_index * 4 + direction];
}

double AbelianGaugeCube::get_field(int nx, int ny, int nz, int nt, int direction) const
{
	return get_field(to_site_index(nx, ny, nz, nt), direction);
}