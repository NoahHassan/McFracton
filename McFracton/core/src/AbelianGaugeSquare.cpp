#include "AbelianGaugeSquare.h"

#include <assert.h>
#include <iostream>

#ifndef PI
#define PI 3.1415926535897932384
#endif

AbelianGaugeSquare::AbelianGaugeSquare(int linear_size, int temporal_size)
	:
	linear_size(linear_size),
	temporal_size(temporal_size),
	nSites(linear_size * linear_size * temporal_size),
	nPlaqs(linear_size * linear_size * temporal_size * 3),
	System(linear_size * linear_size * temporal_size * 3, 0) // A_x, A_y, A_t on each site
{
	site_fields = std::vector<double>(n_site_variables);
	plaq_fields = std::vector<double>(n_site_variables * 3);

	std::random_device rd;
	rng = std::mt19937(rd());
	overrelax_dst = std::uniform_real_distribution<double>(-1.0, 1.0);

	std::for_each(site_fields.begin(), site_fields.end(), [&](double& d) { d = overrelax_dst(rng); });
}

double AbelianGaugeSquare::getEnergy() const
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

double AbelianGaugeSquare::proposeSiteFlip(int index, double angle) const
{
	int site_index = index / 3;
	int type = index % 3;

	int nt = site_index / (linear_size * linear_size);
	int ny = (site_index % (linear_size * linear_size)) / linear_size;
	int nx = site_index % linear_size;

	double current = site_fields[index];

	switch (type)
	{
	case 0:
	{
		return getLocalEnergy_x(nx, ny, nt, current + angle) - getLocalEnergy_x(nx, ny, nt, current);
		break;
	}
	case 1:
	{
		return getLocalEnergy_y(nx, ny, nt, current + angle) - getLocalEnergy_y(nx, ny, nt, current);
		break;
	}
	case 2:
	{
		return getLocalEnergy_t(nx, ny, nt, current + angle) - getLocalEnergy_t(nx, ny, nt, current);
		break;
	}
	default:
	{
		throw("type should not be larger than 2");
		break;
	}
	}

	return 0.0;
}

double AbelianGaugeSquare::proposePlaqFlip(int index, double angle) const
{
	return plaq_fields[index];
}

void AbelianGaugeSquare::UpdateSite(int index, double angle)
{
	site_fields[index] += angle;
}

void AbelianGaugeSquare::UpdatePlaq(int index, double angle)
{
	plaq_fields[index] += angle;
}

void AbelianGaugeSquare::OverrelaxSite(int index)
{
	int site_index = index / 3;
	int direction = site_index % 3;

	int nt = site_index / (linear_size * linear_size);
	int ny = (site_index % (linear_size * linear_size)) / linear_size;
	int nx = site_index % linear_size;

	double random_shift = overrelax_dst(rng);
	// A_i(r) --> A_i(r) + f(r + i) - f(r) = A_i(r) - f(r)
	site_fields[to_site_index(nx, ny, nt) * 3 + 0] += -random_shift;
	site_fields[to_site_index(nx, ny, nt) * 3 + 1] += -random_shift;
	site_fields[to_site_index(nx, ny, nt) * 3 + 2] += -random_shift;

	// A_x-1(r) --> A_x-1(r) + f(r) - f(r-x) = A_x-1(r) + f(r)
	site_fields[to_site_index((nx - 1 + linear_size) % linear_size, ny, nt) * 3 + 0] += random_shift;
	// A_y-1(r) --> A_y-1(r) + f(r) - f(r-y) = A_y-1(r) + f(r)
	site_fields[to_site_index(nx, (ny - 1 + linear_size) % linear_size, nt) * 3 + 1] += random_shift;
	// A_z-1(r) --> A_z-1(r) + f(r) - f(r-z) = A_z-1(r) + f(r)
	site_fields[to_site_index(nx, ny, (nt - 1 + temporal_size) % temporal_size) * 3 + 2] += random_shift;
}

/// <summary>
/// Returns a vector of the integer monople values at each plaquette
/// </summary>
/// <returns></returns>
std::vector<int> AbelianGaugeSquare::getMonopoles() const
{
	std::vector<int> monopoles(nSites);
	for (int n_space = 0; n_space < linear_size * linear_size; n_space++)
	{
		int nx = n_space % linear_size;
		int ny = n_space / linear_size;
		for (int nt = 0; nt < temporal_size; nt++)
		{
			// Because every site corresponds to 3 plaquettes
			std::vector<std::pair<int, int>> cube_faces = {
				{to_site_index(nx, ny, nt), 0},
				{to_site_index(nx, ny, (nt + 1) % temporal_size),0},
				{to_site_index(nx, ny, nt), 1},
				{to_site_index(nx, (ny + 1) % linear_size, nt),1},
				{to_site_index(nx, ny, nt), 2},
				{to_site_index((nx + 1) % linear_size, ny, nt), 2}
			};

			double divergence = 0.0;
			for (int i = 0; i < cube_faces.size(); i++)
			{
				double flux_sign = double(1 - 2 * (i % 2));
				int plaqIndex = cube_faces[i].first * 3 + cube_faces[i].second;
				divergence += mapToCircle(flux_sign * sum_plaquette(getPlaqConnectedFields(plaqIndex)));
			}

			if (divergence >= 1.0 - 1e-5 || divergence <= -1.0 + 1e-5)
			{
				monopoles[to_site_index(nx, ny, nt)] = (int)divergence;
			}
		}
	}
	return monopoles;
}

std::vector<double> AbelianGaugeSquare::getFluxes_z() const
{
	std::vector<double> fluxes(nSites);
	for (int n_space = 0; n_space < linear_size * linear_size; n_space++)
	{
		int nx = n_space % linear_size;
		int ny = n_space / linear_size;
		for (int nt = 0; nt < temporal_size; nt++)
		{
			// Because every site corresponds to 3 plaquettes
			double flux = mapToCircle(
				sum_plaquette(
					getPlaqConnectedFields(
						to_site_index(nx, ny, nt) * 3 + 0
					)
				)
			);

			fluxes[to_site_index(nx, ny, nt)] = flux;
		}
	}

	return fluxes;
}

void AbelianGaugeSquare::LogToFile(std::ofstream& outfile) const
{
	int n_monopoles = 0;
	const auto monopoles = getMonopoles();
	std::for_each(monopoles.begin(), monopoles.end(), [&n_monopoles](int m) { n_monopoles += std::abs(m); });

	outfile << n_monopoles;
}

double AbelianGaugeSquare::getLocalEnergy_x(int nx, int ny, int nt, double angle) const
{
	// Connected to the x-link are four plaquettes: (nx,ny,nt,0),(nx,ny,nt,1),(nx,ny-1,nt,0),(nx,ny,nt-1,1)
	double energy = 0.0;
	std::vector<std::vector<std::pair<int,int>>> connected_plaquettes = {
		getPlaqConnectedFields(nx, ny, nt, 0),
		getPlaqConnectedFields(nx, ny, nt, 1),
		getPlaqConnectedFields(nx, (ny - 1 + linear_size) % linear_size, nt, 0),
		getPlaqConnectedFields(nx, ny, (nt - 1 + temporal_size) % temporal_size, 1)
	};

	for (int n = 0; n < connected_plaquettes.size(); n++)
	{
		auto plaquette = connected_plaquettes[n];
		int angle_index = (4 - n) % 4;
		double plaquette_sum = sum_plaquette(plaquette, angle_index, angle);

		energy += cos(PI * plaquette_sum);
	}

	return -energy;
}

double AbelianGaugeSquare::getLocalEnergy_y(int nx, int ny, int nt, double angle) const
{
	// Connected to the y-link are four plaquettes: (nx,ny,nt,2),(nx-1,ny,nt,0),(nx,ny,nt-1,2),(nx,ny,nt,0)
	double energy = 0.0;
	std::vector<std::vector<std::pair<int, int>>> connected_plaquettes = {
		getPlaqConnectedFields(nx, ny, nt, 2),
		getPlaqConnectedFields((nx - 1 + linear_size) % linear_size, ny, nt, 0),
		getPlaqConnectedFields(nx, ny, (nt - 1 + temporal_size) % temporal_size, 2),
		getPlaqConnectedFields(nx, ny, nt, 0)
	};

	for (int n = 0; n < connected_plaquettes.size(); n++)
	{
		auto plaquette = connected_plaquettes[n];
		int angle_index = n;
		double plaquette_sum = sum_plaquette(plaquette, angle_index, angle);

		energy += cos(PI * plaquette_sum);
	}

	return -energy;
}

double AbelianGaugeSquare::getLocalEnergy_t(int nx, int ny, int nt, double angle) const
{
	// Connected to the t-link are four plaquettes: (nx,ny,nt,1),(nx,ny-1,nt,2),(nx-1,ny,nt,1),(nx,ny,nt,2)
	double energy = 0.0;
	std::vector<std::vector<std::pair<int, int>>> connected_plaquettes = {
		getPlaqConnectedFields(nx, ny, nt, 1),
		getPlaqConnectedFields(nx, (ny - 1 + linear_size) % linear_size, nt, 2),
		getPlaqConnectedFields((nx - 1 + linear_size) % linear_size, ny, nt, 1),
		getPlaqConnectedFields(nx, ny, nt, 2)
	};

	for (int n = 0; n < connected_plaquettes.size(); n++)
	{
		auto plaquette = connected_plaquettes[n];
		int angle_index = n;
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
const std::vector<std::pair<int, int>> AbelianGaugeSquare::getPlaqConnectedFields(int nx, int ny, int nt, int type) const
{
	switch (type)
	{
	case 0: // xy-plaquette with normal in positive t-direction: A_x(r)+A_y(r+x)-A_x(r+y)-A_y(r)
		return {
			{to_site_index(nx, ny, nt), 0},							// A_x(r)
			{to_site_index((nx + 1) % linear_size, ny, nt), 1},		// A_y(r+x)
			{to_site_index(nx, (ny + 1) % linear_size, nt), 0},		// A_x(r+y)
			{to_site_index(nx, ny, nt), 1}							// A_y(r)
		};
		break;
	case 1: // tx-plaquette with normal in positive y-direction: A_t(r)+A_x(r+t)-A_t(r+x)-A_x(r)
		return {
			{to_site_index(nx, ny, nt), 2},							// A_t(r)
			{to_site_index(nx, ny, (nt + 1) % temporal_size), 0},	// A_x(r+t)
			{to_site_index((nx + 1) % linear_size, ny, nt), 2},		// A_t(r+x)
			{to_site_index(nx, ny, nt), 0}							// A_x(r)
		};
		break;
	case 2: // yt-plaquette with normal in positive x-direction: A_y(r)+A_t(r+y)-A_y(r+t)-A_t(r)
		return {
			{to_site_index(nx, ny, nt), 1},							// A_y(r)
			{to_site_index(nx, (ny + 1) % linear_size, nt), 2},		// A_t(r+y)
			{to_site_index(nx, ny, (nt + 1) % temporal_size), 1},	// A_y(r+t)
			{to_site_index(nx, ny, nt), 2}							// A_t(r)
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
const std::vector<std::pair<int, int>> AbelianGaugeSquare::getPlaqConnectedFields(int plaqIndex) const
{
	// Each site connects uniquely to three plaquettes
	// Hence plaqIndex / 3 is the site and plaqIndex % 3 is the type
	int site_index = plaqIndex / 3;
	int plaq_type = plaqIndex % 3;

	int nt = site_index / (linear_size * linear_size);
	int ny = (site_index % (linear_size * linear_size)) / linear_size;
	int nx = site_index % linear_size;

	return getPlaqConnectedFields(nx, ny, nt, plaq_type);
}

// e.g. A_x(r) + A_y(r+x) - A_x(r+y) - A_y(r)
const double AbelianGaugeSquare::sum_plaquette(const std::vector<std::pair<int, int>>& plaquette) const
{
	double sum = 0.0;
	for (int i = 0; i < plaquette.size(); i++)
	{
		double field_sign = double(1 - 2 * (i / 2)); // +1, +1, -1, -1
		sum += field_sign * get_field(plaquette[i].first, plaquette[i].second);
	}
	return sum;
}

const double AbelianGaugeSquare::sum_plaquette(const std::vector<std::pair<int, int>>& plaquette, int angle_index, double angle) const
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

int AbelianGaugeSquare::to_site_index(int nx, int ny, int nt) const
{
	assert(nx < linear_size);
	assert(ny < linear_size);
	assert(nt < temporal_size);
	return (nt * linear_size + ny) * linear_size + nx;
}

double AbelianGaugeSquare::get_field(int site_index, int direction) const
{
	assert(0 <= direction && direction <= 2);
	assert(site_index <= nSites);

	return site_fields[site_index * 3 + direction];
}

double AbelianGaugeSquare::get_field(int nx, int ny, int nt, int direction) const
{
	return get_field(to_site_index(nx, ny, nt), direction);
}

double AbelianGaugeSquare::mapToCircle(const double& d) const
{
	if (d >= 0.5)
		return d - int(d + 0.5);
	else if (d <= 0.5)
		return d - int(d - 0.5);
	else
		return d;
}