#include "AbelianGaugeCube.h"

#include <assert.h>
#include <iostream>

#ifndef PI
#define PI 3.1415926535897932384
#endif

AbelianGaugeCube::AbelianGaugeCube(int linear_size, int temporal_size)
	:
	linear_size(linear_size),
	temporal_size(temporal_size),
	nSites(linear_size * linear_size * linear_size * temporal_size),
	nPlaqs(linear_size * linear_size * linear_size * temporal_size * 6), // xy xz xt yz yt zt
	System(linear_size * linear_size * linear_size * temporal_size * 4, 0) // A_x, A_y, A_z, A_t on each site
{
	site_fields = std::vector<double>(n_site_variables);
	plaq_fields = std::vector<double>(n_site_variables * 3);

	std::random_device rd;
	rng = std::mt19937(rd());
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
		energy += cos(PI * plaquette_sum);
	}

	return -energy;
}

double AbelianGaugeCube::proposeSiteFlip(int index, double angle) const
{
	int site_index = index / 4;
	int type = index % 4;

	const int cube_size = linear_size * linear_size * linear_size;
	const int plane_size = linear_size * linear_size;

	int nt = site_index / cube_size;
	int nz = (site_index % cube_size) / plane_size;
	int ny = (site_index % plane_size) / linear_size;
	int nx = site_index % linear_size;

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
	{
		throw("type should not be larger than 2");
		break;
	}
	}

	return 0.0;
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
	int direction = site_index % 4;

	const int cube_size = linear_size * linear_size * linear_size;
	const int plane_size = linear_size * linear_size;

	int nt = site_index / cube_size;
	int nz = (site_index % cube_size) / plane_size;
	int ny = (site_index % plane_size) / linear_size;
	int nx = site_index % linear_size;

	double random_shift = overrelax_dst(rng);
	// A_i(r) --> A_i(r) + f(r + i) - f(r) = A_i(r) - f(r)
	site_fields[to_site_index(nx, ny, nz, nt) * 4 + 0] += -random_shift;
	site_fields[to_site_index(nx, ny, nz, nt) * 4 + 1] += -random_shift;
	site_fields[to_site_index(nx, ny, nz, nt) * 4 + 2] += -random_shift;
	site_fields[to_site_index(nx, ny, nz, nt) * 4 + 3] += -random_shift;

	// A_x-1(r) --> A_x-1(r) + f(r) - f(r-x) = A_x-1(r) + f(r)
	site_fields[to_site_index((nx - 1 + linear_size) % linear_size, ny, nz, nt) * 4 + 0] += random_shift;
	// A_y-1(r) --> A_y-1(r) + f(r) - f(r-y) = A_y-1(r) + f(r)
	site_fields[to_site_index(nx, (ny - 1 + linear_size) % linear_size, nz, nt) * 4 + 1] += random_shift;
	// A_z-1(r) --> A_z-1(r) + f(r) - f(r-z) = A_z-1(r) + f(r)
	site_fields[to_site_index(nx, ny, (nz - 1 + linear_size) % linear_size, nt) * 4 + 2] += random_shift;
	// A_t-1(r) --> A_t-1(r) + f(r) - f(r-t) = A_t-1(r) + f(r)
	site_fields[to_site_index(nx, ny, nz, (nt - 1 + temporal_size) % temporal_size) * 4 + 3] += random_shift;
}

/// <summary>
/// Returns a vector of the integer monople values at each plaquette
/// </summary>
/// <returns></returns>
std::vector<int> AbelianGaugeCube::getMonopoles() const
{
	std::vector<int> monopoles(nSites*8);
	const int cube_size = linear_size * linear_size * linear_size;
	const int plane_size = linear_size * linear_size;
	for (int n_space = 0; n_space < cube_size; n_space++)
	{
		for (int nt = 0; nt < temporal_size; nt++)
		{
			int nz = n_space / plane_size;
			int ny = (n_space % plane_size) / linear_size;
			int nx = n_space % linear_size;

			// Because every site corresponds to 6 plaquettes
			// 0:xy, 1:xz, 2:xt, 3:yz, 4:yt, 5:zt
			std::vector<std::pair<int, int>> hypercube_faces = {
				// connected to (0,0,0,0) [all positive]
				{to_site_index(nx, ny, nz, nt), 0}, // xy 0
				{to_site_index(nx, ny, nz, nt), 1}, // xz 1
				{to_site_index(nx, ny, nz, nt), 2}, // xt 2
				{to_site_index(nx, ny, nz, nt), 3}, // yz 3
				{to_site_index(nx, ny, nz, nt), 4}, // yt 4
				{to_site_index(nx, ny, nz, nt), 5}, // zt 5

				// [all positive]
				// connected to (1,1,0,0)
				{to_site_index((nx + 1) % linear_size, (ny + 1) % linear_size, nz, nt), 5},   // zt 6

				// connected to (1,0,1,0)
				{to_site_index((nx + 1) % linear_size, ny, (nz + 1) % linear_size, nt), 4},   // yt 7

				// connected to (1,0,0,1)
				{to_site_index((nx + 1) % linear_size, ny, nz, (nt + 1) % temporal_size), 3}, // yz 8

				// connected to (0,1,1,0)
				{to_site_index(nx, (ny + 1) % linear_size, (nz + 1) % linear_size, nt), 2},	  // xt 9

				// connected to (0,1,0,1)
				{to_site_index(nx, (ny + 1) % linear_size, nz, (nt + 1) % temporal_size), 1}, // xz 10

				// connected to (0,0,1,1)
				{to_site_index(nx, ny, (nz + 1) % linear_size, (nt + 1) % temporal_size), 0}, // xy 11

				// [all negative]
				// connected to (1,0,0,0)
				{to_site_index((nx + 1) % linear_size, ny, nz, nt), 3}, // yz 12
				{to_site_index((nx + 1) % linear_size, ny, nz, nt), 4}, // yt 13
				{to_site_index((nx + 1) % linear_size, ny, nz, nt), 5}, // zt 14

				// connected to (0,1,0,0)
				{to_site_index(nx, (ny + 1) % linear_size, nz, nt), 1}, // xz 15
				{to_site_index(nx, (ny + 1) % linear_size, nz, nt), 2}, // xt 16
				{to_site_index(nx, (ny + 1) % linear_size, nz, nt), 5}, // zt 17

				// connected to (0,0,1,0)
				{to_site_index(nx, ny, (nz + 1) % linear_size, nt), 0}, // xy 18
				{to_site_index(nx, ny, (nz + 1) % linear_size, nt), 2}, // xt 19
				{to_site_index(nx, ny, (nz + 1) % linear_size, nt), 4}, // yt 20

				// connected to (0,0,0,1)
				{to_site_index(nx, ny, nz, (nt + 1) % temporal_size), 0}, // xy 21
				{to_site_index(nx, ny, nz, (nt + 1) % temporal_size), 1}, // xz 22
				{to_site_index(nx, ny, nz, (nt + 1) % temporal_size), 3}  // yz 23
			};

			std::vector<std::vector<int>> cubes = {
				{0,1,3,18,15,12},	// xy xz yz, xy xz yz (+ - +)
				{21,22,23,11,10,8}, // xy xz yz, xy xz yz (+ - +)

				{2,1,5,19,22,14},	// xt xz zt, xt xz zt (+ - +)
				{16,15,17,9,10,6},	// xt xz zt, xt xz zt (+ - +)
				{4,3,5,20,23,17},	// yt yz zt, yt yz zt (+ - +)
				{13,8,14,7,12,6},	// yt yz zt, yt yz zt (+ - +)
				{0,1,4,21,16,13},	// xy xt yt, xy xt yt (+ - +)
				{11,19,20,18,9,7}	// xy xt yt, xy xt yt (+ - +)
			};

			// loop through all cubes of the hypercube boundary
			for (int i = 0; i < 8; i++)
			{
				double divergence = 0.0;

				// get individual cube of the boundary
				std::vector<int> cube_face_indices = cubes[i];
				
				// loop over all faces of the cube
				for (int n_face = 0; n_face < cube_face_indices.size(); n_face++)
				{
					// flux sign depends on front/back, left/right, top/bottom
					// according to the curl definitions below and the ordering in the
					// struct above the signs must be + - + - + -
					double flux_sign = double(1 - 2 * (i % 2));

					// get actual face in spacetime
					auto hypercube_face = hypercube_faces[cube_face_indices[n_face]];
					int plaqIndex = hypercube_face.first * 6 + hypercube_face.second;

					// calculate the flux
					divergence += mapToCircle(flux_sign * sum_plaquette(getPlaqConnectedFields(plaqIndex)));
				}

				if (divergence >= 1.0 - 1e-5 || divergence <= -1.0 + 1e-5)
				{
					monopoles[to_site_index(nx, ny, nz, nt) + i] = (int)divergence;
				}
			}
		}
	}
	return monopoles;
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
			double flux = mapToCircle(
				sum_plaquette(
					getPlaqConnectedFields(
						to_site_index(nx, ny, nz, nt) * 4 + 0
					)
				)
			);

			fluxes[to_site_index(nx, ny, nz, nt)] = flux;
		}
	}

	return fluxes;
}

//void AbelianGaugeSquare::LogToFile(std::ofstream& outfile) const
//{
//	//int n_monopoles = 0;
//	//const auto monopoles = getMonopoles();
//	//std::for_each(monopoles.begin(), monopoles.end(), [&n_monopoles](int m) { n_monopoles += std::abs(m); });
//
//	outfile << getEnergy();
//}

System::Observables AbelianGaugeCube::Measure(double T) const
{
	System::Observables observables;
	observables.energy = getEnergy();
	observables.helicity_modulus = 0.0;

	int n_a = 0;
	int n_b = 0;
	const auto monopoles = getMonopoles();
	for (int n = 0; n < monopoles.size(); n++)
	{
		if (monopoles[n] < 0)
			n_b++;
		else if (monopoles[n] > 0)
			n_a++;
	}

	observables.n_defects_a = n_a;
	observables.n_defects_b = n_b;

	int plane_size = linear_size * linear_size;
	int cube_size = plane_size * linear_size;
	std::vector<double> loops(cube_size);
	for (int n_space = 0; n_space < cube_size; n_space++)
	{
		double field_sum = 0.0;
		for (int nt = 0; nt < temporal_size; nt++)
		{
			int site_index = cube_size * nt + n_space;
			field_sum += get_field(site_index, 4);
		}
		loops[n_space] = cos(field_sum);
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
		getPlaqConnectedFields(nx, (ny - 1 + linear_size) % linear_size, nz, nt, 0),	// xy
		getPlaqConnectedFields(nx, ny, (nz - 1 + linear_size) % linear_size, nt, 1),	// xz
		getPlaqConnectedFields(nx, ny, nz, (nt - 1 + temporal_size) % temporal_size, 2) // xt
	};

	for (int n = 0; n < connected_plaquettes.size(); n++)
	{
		auto plaquette = connected_plaquettes[n];
		int angle_index = (n / 3) * 2; // because the links appear at first or third position
		double plaquette_sum = sum_plaquette(plaquette, angle_index, angle);

		energy += cos(PI * plaquette_sum);
	}

	return -energy;
}

double AbelianGaugeCube::getLocalEnergy_y(int nx, int ny, int nz, int nt, double angle) const
{
	// 0: xy, 1: xz, 2: xt, 3: yz, 4: yt, 5: zt
	double energy = 0.0;
	std::vector<std::vector<std::pair<int, int>>> connected_plaquettes = {
		getPlaqConnectedFields(nx, ny, nz, nt, 0),										// xy
		getPlaqConnectedFields((nx - 1 + linear_size) % linear_size, ny, nz, nt, 0),	// xy
		getPlaqConnectedFields(nx, ny, nz, nt, 3),										// yz
		getPlaqConnectedFields(nx, ny, (nz - 1 + linear_size) % linear_size, nt, 3),	// yz
		getPlaqConnectedFields(nx, ny, nz, nt, 4),										// yt
		getPlaqConnectedFields(nx, ny, nz, (nt - 1 + temporal_size) % temporal_size, 4)	// yt
	};

	int angle_indices[] = { 3, 1, 0, 2, 0, 2 };
	for (int n = 0; n < connected_plaquettes.size(); n++)
	{
		auto plaquette = connected_plaquettes[n];
		int angle_index = angle_indices[n];
		double plaquette_sum = sum_plaquette(plaquette, angle_index, angle);

		energy += cos(PI * plaquette_sum);
	}

	return -energy;
}

double AbelianGaugeCube::getLocalEnergy_z(int nx, int ny, int nz, int nt, double angle) const
{
	// 0: xy, 1: xz, 2: xt, 3: yz, 4: yt, 5: zt
	double energy = 0.0;
	std::vector<std::vector<std::pair<int, int>>> connected_plaquettes = {
		getPlaqConnectedFields(nx, ny, nz, nt, 1),										// xz
		getPlaqConnectedFields((nx - 1 + linear_size) % linear_size, ny, nz, nt, 1),	// xz
		getPlaqConnectedFields(nx, ny, nz, nt, 3),										// yz
		getPlaqConnectedFields(nx, (ny - 1 + linear_size) % linear_size, nz, nt, 3),	// yz
		getPlaqConnectedFields(nx, ny, nz, nt, 5),										// zt
		getPlaqConnectedFields(nx, ny, nz, (nt - 1 + temporal_size) % temporal_size, 5)	// zt
	};

	int angle_indices[] = { 3, 1, 3, 1, 0, 2 };
	for (int n = 0; n < connected_plaquettes.size(); n++)
	{
		auto plaquette = connected_plaquettes[n];
		int angle_index = angle_indices[n];
		double plaquette_sum = sum_plaquette(plaquette, angle_index, angle);

		energy += cos(PI * plaquette_sum);
	}

	return -energy;
}

double AbelianGaugeCube::getLocalEnergy_t(int nx, int ny, int nz, int nt, double angle) const
{
	// 0: xy, 1: xz, 2: xt, 3: yz, 4: yt, 5: zt
	double energy = 0.0;
	std::vector<std::vector<std::pair<int, int>>> connected_plaquettes = {
		getPlaqConnectedFields(nx, ny, nz, nt, 2),										// xt
		getPlaqConnectedFields((nx - 1 + linear_size) % linear_size, ny, nz, nt, 2),	// xt
		getPlaqConnectedFields(nx, ny, nz, nt, 4),										// yt
		getPlaqConnectedFields(nx, (ny - 1 + linear_size) % linear_size, nz, nt, 4),	// yt
		getPlaqConnectedFields(nx, ny, nz, nt, 5),										// zt
		getPlaqConnectedFields(nx, ny, (nz - 1 + linear_size) % linear_size, nt, 5)		// zt
	};

	int angle_indices[] = { 3, 1, 3, 1, 3, 1 };
	for (int n = 0; n < connected_plaquettes.size(); n++)
	{
		auto plaquette = connected_plaquettes[n];
		int angle_index = angle_indices[n];
		double plaquette_sum = sum_plaquette(plaquette, angle_index, angle);

		energy += cos(PI * plaquette_sum);
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
const std::vector<std::pair<int, int>> AbelianGaugeCube::getPlaqConnectedFields(int nx, int ny, int nz, int nt, int type) const
{
	// 0:xy, 1:xz, 2:xt, 3:yz, 4:yt, 5:zt
	switch (type)
	{
	case 0: // xy-plaquette: A_x(r)+A_y(r+x)-A_x(r+y)-A_y(r)
		return {
			{to_site_index(nx, ny, nz, nt), 0},							// A_x(r)
			{to_site_index((nx + 1) % linear_size, ny, nz, nt), 1},		// A_y(r+x)
			{to_site_index(nx, (ny + 1) % linear_size, nz, nt), 0},		// A_x(r+y)
			{to_site_index(nx, ny, nz, nt), 1}							// A_y(r)
		};
		break;
	case 1: // xz-plaquette: A_x(r)+A_z(r+x)-A_x(r+z)-A_z(r)
		return {
			{to_site_index(nx, ny, nz, nt), 0},							// A_x(r)
			{to_site_index((nx + 1) % linear_size, ny, nz, nt), 2},		// A_z(r+x)
			{to_site_index(nx, ny, (nz + 1) % linear_size, nt), 0},		// A_x(r+z)
			{to_site_index(nx, ny, nz, nt), 2}							// A_z(r)
		};
		break;
	case 2: // xt-plaquette: A_x(r)+A_t(r+x)-A_x(r+t)-A_t(r)
		return {
			{to_site_index(nx, ny, nz, nt), 0},							// A_x(r)
			{to_site_index((nx + 1) % linear_size, ny, nz, nt), 3},		// A_t(r+x)
			{to_site_index(nx, ny, nz, (nt + 1) % temporal_size), 0},	// A_x(r+t)
			{to_site_index(nx, ny, nz, nt), 3}							// A_t(r)
		};
		break;
	case 3: // yz-plaquette: A_y(r)+A_z(r+y)-A_y(r+z)-A_z(r)
		return {
			{to_site_index(nx, ny, nz, nt), 1},							// A_y(r)
			{to_site_index(nx, (ny + 1) % linear_size, nz, nt), 2},		// A_z(r+y)
			{to_site_index(nx, ny, (nz + 1) % linear_size, nt), 1},		// A_x(r+z)
			{to_site_index(nx, ny, nz, nt), 2}							// A_z(r)
		};
		break;
	case 4: // yt-plaquette: A_y(r)+A_t(r+y)-A_y(r+t)-A_t(r)
		return {
			{to_site_index(nx, ny, nz, nt), 1},							// A_y(r)
			{to_site_index(nx, (ny + 1) % linear_size, nz, nt), 3},		// A_t(r+y)
			{to_site_index(nx, ny, nz, (nt + 1) % temporal_size), 1},	// A_y(r+t)
			{to_site_index(nx, ny, nz, nt), 3}							// A_t(r)
		};
		break;
	case 5: // zt-plaquette: A_z(r)+A_t(r+z)-A_z(r+t)-A_t(r)
		return {
			{to_site_index(nx, ny, nz, nt), 2},							// A_z(r)
			{to_site_index(nx, ny, (nz + 1) % linear_size, nt), 3},		// A_t(r+z)
			{to_site_index(nx, ny, nz, (nt + 1) % temporal_size), 2},	// A_z(r+t)
			{to_site_index(nx, ny, nz, nt), 3}							// A_t(r)
		};
		break;
	default:
		throw("No plaquette this type defined");
		break;
	}
}

/// <summary>
/// Given a plaquette, returns the corresponding fields as pairs of site and direction.
/// </summary>
/// <param name="plaqIndex"></param>
/// <returns></returns>
const std::vector<std::pair<int, int>> AbelianGaugeCube::getPlaqConnectedFields(int plaqIndex) const
{
	// Each site connects uniquely to three plaquettes
	// Hence plaqIndex / 3 is the site and plaqIndex % 3 is the type
	int site_index = plaqIndex / 6;
	int plaq_type = plaqIndex % 6;

	const int cube_size = linear_size * linear_size * linear_size;
	const int plane_size = linear_size * linear_size;

	int nt = site_index / cube_size;
	int nz = (site_index % cube_size) / plane_size;
	int ny = (site_index % plane_size) / linear_size;
	int nx = site_index % linear_size;

	return getPlaqConnectedFields(nx, ny, nz, nt, plaq_type);
}

// e.g. A_x(r) + A_y(r+x) - A_x(r+y) - A_y(r)
const double AbelianGaugeCube::sum_plaquette(const std::vector<std::pair<int, int>>& plaquette) const
{
	double sum = 0.0;
	for (int i = 0; i < plaquette.size(); i++)
	{
		double field_sign = double(1 - 2 * (i / 2)); // +1, +1, -1, -1
		sum += field_sign * get_field(plaquette[i].first, plaquette[i].second);
	}
	return sum;
}

const double AbelianGaugeCube::sum_plaquette(const std::vector<std::pair<int, int>>& plaquette, int angle_index, double angle) const
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
	assert(nx < linear_size);
	assert(ny < linear_size);
	assert(nz < linear_size);
	assert(nt < temporal_size);
	return ((nt * linear_size + nz) * linear_size + ny) * linear_size + nx;
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

double AbelianGaugeCube::mapToCircle(const double& d) const
{
	if (d >= 0.5)
		return d - int(d + 0.5);
	else if (d <= 0.5)
		return d - int(d - 0.5);
	else
		return d;
}