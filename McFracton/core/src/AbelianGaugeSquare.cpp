#include "AbelianGaugeSquare.h"

#include <assert.h>
#include <iostream>
#include <stdexcept>

#include "MathUtil.h"

using mcf::kPi;

AbelianGaugeSquare::AbelianGaugeSquare(int linear_size, int temporal_size, unsigned int seed)
	:
	// A_x, A_y, A_t on each site
	System(mcf::PeriodicLattice({ linear_size, linear_size, temporal_size }, { "x", "y", "t" }),
		linear_size * linear_size * temporal_size * 3),
	linear_size(linear_size),
	temporal_size(temporal_size),
	nSites(linear_size * linear_size * temporal_size),
	nPlaqs(linear_size * linear_size * temporal_size * 3)
{
	// seed == 0 means "pick a fresh, unpredictable seed"; any other value is reproducible.
	rng = std::mt19937(seed != 0 ? seed : std::random_device{}());
	overrelax_dst = std::uniform_real_distribution<double>(-1.0, 1.0);

	std::for_each(fields.begin(), fields.end(), [&](double& d) { d = overrelax_dst(rng); });
}

double AbelianGaugeSquare::getEnergy() const
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

double AbelianGaugeSquare::proposeUpdate(int index, double delta) const
{
	int site_index = index / 3;
	int type = index % 3;

	int nt = lattice.coord(site_index, 2);
	int ny = lattice.coord(site_index, 1);
	int nx = lattice.coord(site_index, 0);

	double current = fields[index];

	switch (type)
	{
	case 0:
	{
		return getLocalEnergy_x(nx, ny, nt, current + delta) - getLocalEnergy_x(nx, ny, nt, current);
		break;
	}
	case 1:
	{
		return getLocalEnergy_y(nx, ny, nt, current + delta) - getLocalEnergy_y(nx, ny, nt, current);
		break;
	}
	case 2:
	{
		return getLocalEnergy_t(nx, ny, nt, current + delta) - getLocalEnergy_t(nx, ny, nt, current);
		break;
	}
	default:
		throw std::logic_error("AbelianGaugeSquare: link direction must be 0, 1 or 2");
	}
}

void AbelianGaugeSquare::overrelax(int index)
{
	int site_index = index / 3;

	int nt = lattice.coord(site_index, 2);
	int ny = lattice.coord(site_index, 1);
	int nx = lattice.coord(site_index, 0);

	double random_shift = overrelax_dst(rng);
	// A_i(r) --> A_i(r) + f(r + i) - f(r) = A_i(r) - f(r)
	fields[to_site_index(nx, ny, nt) * 3 + 0] += -random_shift;
	fields[to_site_index(nx, ny, nt) * 3 + 1] += -random_shift;
	fields[to_site_index(nx, ny, nt) * 3 + 2] += -random_shift;

	// A_x-1(r) --> A_x-1(r) + f(r) - f(r-x) = A_x-1(r) + f(r)
	fields[to_site_index(lattice.wrap(0, nx - 1), ny, nt) * 3 + 0] += random_shift;
	// A_y-1(r) --> A_y-1(r) + f(r) - f(r-y) = A_y-1(r) + f(r)
	fields[to_site_index(nx, lattice.wrap(1, ny - 1), nt) * 3 + 1] += random_shift;
	// A_z-1(r) --> A_z-1(r) + f(r) - f(r-z) = A_z-1(r) + f(r)
	fields[to_site_index(nx, ny, lattice.wrap(2, nt - 1)) * 3 + 2] += random_shift;
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
				{to_site_index(nx, ny, lattice.wrap(2, nt + 1)),0},
				{to_site_index(nx, ny, nt), 1},
				{to_site_index(nx, lattice.wrap(1, ny + 1), nt),1},
				{to_site_index(nx, ny, nt), 2},
				{to_site_index(lattice.wrap(0, nx + 1), ny, nt), 2}
			};

			double divergence = 0.0;
			for (int i = 0; i < cube_faces.size(); i++)
			{
				double flux_sign = double(1 - 2 * (i % 2));
				int plaqIndex = cube_faces[i].first * 3 + cube_faces[i].second;
				divergence += mcf::mapToCircle(flux_sign * sum_plaquette(getPlaqConnectedFields(plaqIndex)));
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
			double flux = mcf::mapToCircle(
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

std::vector<std::string> AbelianGaugeSquare::observableNames() const
{
	return { "Energy", "defects_a", "defects_b" };
}

std::vector<double> AbelianGaugeSquare::measure(double temperature) const
{
	const double energy = getEnergy();

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

	return { energy, (double)n_a, (double)n_b };
}

std::vector<Channel> AbelianGaugeSquare::channels() const
{
	return { { "A_x", ChannelKind::Angle },
			 { "A_y", ChannelKind::Angle },
			 { "A_t", ChannelKind::Angle },
			 { "flux_xy", ChannelKind::Signed },
			 { "monopoles", ChannelKind::Integer } };
}

void AbelianGaugeSquare::fillChannel(int channel, std::vector<double>& out) const
{
	out.assign(lattice.size(), 0.0);

	if (0 <= channel && channel <= 2)
	{
		for (int site = 0; site < lattice.size(); site++)
			out[site] = get_field(site, channel);
		return;
	}
	if (channel == 3)
	{
		out = getFluxes_z();
		return;
	}

	const std::vector<int> monopoles = getMonopoles();
	for (int site = 0; site < lattice.size(); site++)
		out[site] = (double)monopoles[site];
}

double AbelianGaugeSquare::getLocalEnergy_x(int nx, int ny, int nt, double angle) const
{
	// Connected to the x-link are four plaquettes: (nx,ny,nt,0),(nx,ny,nt,1),(nx,ny-1,nt,0),(nx,ny,nt-1,1)
	double energy = 0.0;
	std::vector<std::vector<std::pair<int,int>>> connected_plaquettes = {
		getPlaqConnectedFields(nx, ny, nt, 0),
		getPlaqConnectedFields(nx, ny, nt, 1),
		getPlaqConnectedFields(nx, lattice.wrap(1, ny - 1), nt, 0),
		getPlaqConnectedFields(nx, ny, lattice.wrap(2, nt - 1), 1)
	};

	for (int n = 0; n < connected_plaquettes.size(); n++)
	{
		auto plaquette = connected_plaquettes[n];
		int angle_index = (4 - n) % 4;
		double plaquette_sum = sum_plaquette(plaquette, angle_index, angle);

		energy += cos(kPi * plaquette_sum);
	}

	return -energy;
}

double AbelianGaugeSquare::getLocalEnergy_y(int nx, int ny, int nt, double angle) const
{
	// Connected to the y-link are four plaquettes: (nx,ny,nt,2),(nx-1,ny,nt,0),(nx,ny,nt-1,2),(nx,ny,nt,0)
	double energy = 0.0;
	std::vector<std::vector<std::pair<int, int>>> connected_plaquettes = {
		getPlaqConnectedFields(nx, ny, nt, 2),
		getPlaqConnectedFields(lattice.wrap(0, nx - 1), ny, nt, 0),
		getPlaqConnectedFields(nx, ny, lattice.wrap(2, nt - 1), 2),
		getPlaqConnectedFields(nx, ny, nt, 0)
	};

	for (int n = 0; n < connected_plaquettes.size(); n++)
	{
		auto plaquette = connected_plaquettes[n];
		int angle_index = n;
		double plaquette_sum = sum_plaquette(plaquette, angle_index, angle);

		energy += cos(kPi * plaquette_sum);
	}

	return -energy;
}

double AbelianGaugeSquare::getLocalEnergy_t(int nx, int ny, int nt, double angle) const
{
	// Connected to the t-link are four plaquettes: (nx,ny,nt,1),(nx,ny-1,nt,2),(nx-1,ny,nt,1),(nx,ny,nt,2)
	double energy = 0.0;
	std::vector<std::vector<std::pair<int, int>>> connected_plaquettes = {
		getPlaqConnectedFields(nx, ny, nt, 1),
		getPlaqConnectedFields(nx, lattice.wrap(1, ny - 1), nt, 2),
		getPlaqConnectedFields(lattice.wrap(0, nx - 1), ny, nt, 1),
		getPlaqConnectedFields(nx, ny, nt, 2)
	};

	for (int n = 0; n < connected_plaquettes.size(); n++)
	{
		auto plaquette = connected_plaquettes[n];
		int angle_index = n;
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
std::vector<std::pair<int, int>> AbelianGaugeSquare::getPlaqConnectedFields(int nx, int ny, int nt, int type) const
{
	switch (type)
	{
	case 0: // xy-plaquette with normal in positive t-direction: A_x(r)+A_y(r+x)-A_x(r+y)-A_y(r)
		return {
			{to_site_index(nx, ny, nt), 0},							// A_x(r)
			{to_site_index(lattice.wrap(0, nx + 1), ny, nt), 1},		// A_y(r+x)
			{to_site_index(nx, lattice.wrap(1, ny + 1), nt), 0},		// A_x(r+y)
			{to_site_index(nx, ny, nt), 1}							// A_y(r)
		};
		break;
	case 1: // tx-plaquette with normal in positive y-direction: A_t(r)+A_x(r+t)-A_t(r+x)-A_x(r)
		return {
			{to_site_index(nx, ny, nt), 2},							// A_t(r)
			{to_site_index(nx, ny, lattice.wrap(2, nt + 1)), 0},	// A_x(r+t)
			{to_site_index(lattice.wrap(0, nx + 1), ny, nt), 2},		// A_t(r+x)
			{to_site_index(nx, ny, nt), 0}							// A_x(r)
		};
		break;
	case 2: // yt-plaquette with normal in positive x-direction: A_y(r)+A_t(r+y)-A_y(r+t)-A_t(r)
		return {
			{to_site_index(nx, ny, nt), 1},							// A_y(r)
			{to_site_index(nx, lattice.wrap(1, ny + 1), nt), 2},		// A_t(r+y)
			{to_site_index(nx, ny, lattice.wrap(2, nt + 1)), 1},	// A_y(r+t)
			{to_site_index(nx, ny, nt), 2}							// A_t(r)
		};
		break;
	default:
		throw std::logic_error("AbelianGaugeSquare: plaquette type must be 0, 1 or 2");
	}
}

/// <summary>
/// Given a plaquette, returns the corresponding fields as pairs of site and direction.
/// </summary>
/// <param name="plaqIndex"></param>
/// <returns></returns>
std::vector<std::pair<int, int>> AbelianGaugeSquare::getPlaqConnectedFields(int plaqIndex) const
{
	// Each site connects uniquely to three plaquettes
	// Hence plaqIndex / 3 is the site and plaqIndex % 3 is the type
	int site_index = plaqIndex / 3;
	int plaq_type = plaqIndex % 3;

	return getPlaqConnectedFields(lattice.coord(site_index, 0),
		lattice.coord(site_index, 1),
		lattice.coord(site_index, 2),
		plaq_type);
}

// e.g. A_x(r) + A_y(r+x) - A_x(r+y) - A_y(r)
double AbelianGaugeSquare::sum_plaquette(const std::vector<std::pair<int, int>>& plaquette) const
{
	double sum = 0.0;
	for (int i = 0; i < plaquette.size(); i++)
	{
		double field_sign = double(1 - 2 * (i / 2)); // +1, +1, -1, -1
		sum += field_sign * get_field(plaquette[i].first, plaquette[i].second);
	}
	return sum;
}

double AbelianGaugeSquare::sum_plaquette(const std::vector<std::pair<int, int>>& plaquette, int angle_index, double angle) const
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
	return lattice.index(nx, ny, nt);
}

double AbelianGaugeSquare::get_field(int site_index, int direction) const
{
	assert(0 <= direction && direction <= 2);
	assert(site_index <= nSites);

	return fields[site_index * 3 + direction];
}

double AbelianGaugeSquare::get_field(int nx, int ny, int nt, int direction) const
{
	return get_field(to_site_index(nx, ny, nt), direction);
}