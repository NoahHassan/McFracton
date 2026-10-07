#include "Spiderweb_corrected.h"

#include <assert.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "MathUtil.h"

using mcf::kPi;

namespace {

	// The checkerboard only closes periodically on an even linear size, so anything else is refused
	// here, before the base class sizes its fields. The smallest lattices are refused as well: at
	// linear_size 2 (r+x == r-x) and temporal_size 1 (r-t == r) two neighbours of a site are the
	// same site, and proposeUpdate relies on them all being different.
	int checkedLinearSize(int linear_size, int temporal_size)
	{
		if (linear_size < 4 || linear_size % 2 != 0)
			throw std::invalid_argument("Spiderweb_corrected needs an even linear_size >= 4 (two-site basis)");
		if (temporal_size < 2)
			throw std::invalid_argument("Spiderweb_corrected needs temporal_size >= 2");
		return linear_size;
	}

}

Spiderweb_corrected::Spiderweb_corrected(int linear_size, int temporal_size, double KU)
	:
	System(mcf::PeriodicLattice({ checkedLinearSize(linear_size, temporal_size), linear_size, temporal_size }, { "x", "y", "t" }),
		linear_size * linear_size * temporal_size * 3 / 2), // Because each unit-cell has 3 fields and linear_size counts both crossed and empty squares
	linear_size(linear_size),
	temporal_size(temporal_size),
	nSites(linear_size * linear_size * temporal_size),
	KU(KU)
{
	if (!(KU > 0.0))
		throw std::invalid_argument("Spiderweb_corrected needs KU > 0");
}

void Spiderweb_corrected::randomize(std::mt19937& rng)
{
	std::uniform_real_distribution<double> dst(-1.0, 1.0);
	std::for_each(fields.begin(), fields.end(), [&](double& d) { d = dst(rng); });
}

void Spiderweb_corrected::setGroundState()
{
	std::fill(fields.begin(), fields.end(), 0.0);
}

bool Spiderweb_corrected::isType1(int site_index) const
{
	return (lattice.coord(site_index, 0) + lattice.coord(site_index, 1)) % 2 == 0;
}

// To avoid having two fields on type 1 (A^xy and A^0) and one field on type 2 (A^xx)
// A layout is chosen by which first all spatial components are listed and the
// A0 fields come after in the array
int Spiderweb_corrected::a0Index(int site_index) const
{
	assert(isType1(site_index)); // Important because A^0 lives on links connecting crossed sties
	return nSites + site_index / 2;
}

int Spiderweb_corrected::shifted(int site_index, int dx, int dy, int dt) const
{
	int site = site_index;
	if (dx != 0) site = lattice.neighbor(site, 0, dx);
	if (dy != 0) site = lattice.neighbor(site, 1, dy);
	if (dt != 0) site = lattice.neighbor(site, 2, dt);
	return site;
}

// B(r) = A^xy(r+x) - A^xx(r+x+y) - A^xy(r+y) + A^xx(r-x+y)
//      + A^xy(r-x) - A^xx(r-x-y) - A^xy(r-y) + A^xx(r+x-y)
Spiderweb_corrected::Stencil Spiderweb_corrected::magneticStencil(int r) const
{
	assert(!isType1(r));
	Stencil s;

	s.add(tensorIndex(shifted(r, +1, 0, 0)), +1.0);
	s.add(tensorIndex(shifted(r, +1, +1, 0)), -1.0);
	s.add(tensorIndex(shifted(r, 0, +1, 0)), -1.0);
	s.add(tensorIndex(shifted(r, -1, +1, 0)), +1.0);
	s.add(tensorIndex(shifted(r, -1, 0, 0)), +1.0);
	s.add(tensorIndex(shifted(r, -1, -1, 0)), -1.0);
	s.add(tensorIndex(shifted(r, 0, -1, 0)), -1.0);
	s.add(tensorIndex(shifted(r, +1, -1, 0)), +1.0);
	return s;
}

// E^xy(r) = A^0(r+x+y) - A^0(r-x+y) + A^0(r-x-y) - A^0(r+x-y) - A^xy(r+t) + A^xy(r)
Spiderweb_corrected::Stencil Spiderweb_corrected::electricStencil_xy(int r) const
{
	assert(isType1(r)); // E^xy lives on crossed (type 1) squares.
	Stencil s;

	s.add(a0Index(shifted(r, +1, +1, 0)), +1.0);
	s.add(a0Index(shifted(r, -1, +1, 0)), -1.0);
	s.add(a0Index(shifted(r, -1, -1, 0)), +1.0);
	s.add(a0Index(shifted(r, +1, -1, 0)), -1.0);
	s.add(tensorIndex(shifted(r, 0, 0, +1)), -1.0);
	s.add(tensorIndex(r), +1.0);
	return s;
}

// E^xx(r) = A^0(r+x) - A^0(r+y) + A^0(r-x) - A^0(r-y) - A^xx(r+t) + A^xx(r)
Spiderweb_corrected::Stencil Spiderweb_corrected::electricStencil_xx(int r) const
{
	assert(!isType1(r)); // E^xx lives on empty (type 2) squares.
	Stencil s;

	s.add(a0Index(shifted(r, +1, 0, 0)), +1.0);
	s.add(a0Index(shifted(r, 0, +1, 0)), -1.0);
	s.add(a0Index(shifted(r, -1, 0, 0)), +1.0);
	s.add(a0Index(shifted(r, 0, -1, 0)), -1.0);
	s.add(tensorIndex(shifted(r, 0, 0, +1)), -1.0);
	s.add(tensorIndex(r), +1.0);
	return s;
}

double Spiderweb_corrected::sum(const Stencil& stencil) const
{
	double total = 0.0;
	for (int i = 0; i < stencil.size; i++)
		total += stencil.coefficient[i] * fields[stencil.index[i]];
	return total;
}

double Spiderweb_corrected::magneticField(int site_index) const
{
	return sum(magneticStencil(site_index));
}

double Spiderweb_corrected::electricField_xx(int site_index) const
{
	return sum(electricStencil_xx(site_index));
}

double Spiderweb_corrected::electricField_xy(int site_index) const
{
	return sum(electricStencil_xy(site_index));
}

std::vector<double> Spiderweb_corrected::getLocalEnergies() const
{
	std::vector<double> localEnergies(nSites);
	for (int site = 0; site < nSites; site++)
	{
		if (isType1(site))
			localEnergies[site] = -cos(kPi * electricField_xy(site)) / KU;
		else
			localEnergies[site] = -cos(kPi * magneticField(site)) - cos(kPi * electricField_xx(site)) / KU;
	}
	return localEnergies;
}

double Spiderweb_corrected::getEnergy() const
{
	double energy = 0.0;
	for (double local_energy : getLocalEnergies())
		energy += local_energy;
	return energy;
}

std::vector<int> Spiderweb_corrected::getInstantons() const
{
	// Get B and E fields at every site wrapped to the interval [-1,1)
	std::vector<double> b(nSites, 0.0), e(nSites, 0.0);
	for (int site = 0; site < nSites; site++)
	{
		if (isType1(site))
			e[site] = mcf::mapToCircle(electricField_xy(site));
		else
		{
			b[site] = mcf::mapToCircle(magneticField(site));
			e[site] = mcf::mapToCircle(electricField_xx(site));
		}
	}

	// Instantons are computed as violations of the gauss law Dt B + Dxx E^xy - Dxy E^xx = 0
	std::vector<int> instantons(nSites, 0);
	for (int r = 0; r < nSites; r++)
	{
		if (isType1(r))
			continue;

		const double wrapped_gauss = b[shifted(r, 0, 0, +1)] - b[r]
			// Dxx [E^xy]: the four type 1 neighbours
			+ e[shifted(r, +1, 0, 0)] - e[shifted(r, 0, +1, 0)] + e[shifted(r, -1, 0, 0)] - e[shifted(r, 0, -1, 0)]
			// - Dxy [E^xx]: the four diagonal type 2 neighbours
			- (e[shifted(r, +1, +1, 0)] - e[shifted(r, -1, +1, 0)] + e[shifted(r, -1, -1, 0)] - e[shifted(r, +1, -1, 0)]);
		instantons[r] = (int)std::lround(wrapped_gauss / 2.0);
	}
	return instantons;
}

double Spiderweb_corrected::proposeUpdate(int index, double delta) const
{
	double d_energy = 0.0;

	// Given a stencil (magnetic, electric_xx, electric_xy), checks whether the field at index
	// is inside that stencil and adds the corresponding energy difference to d_energy
	auto accumulate = [&](const Stencil& stencil, double prefactor)
		{
			double old_sum = 0.0;
			double weight = 0.0;
			for (int i = 0; i < stencil.size; i++)
			{
				old_sum += stencil.coefficient[i] * fields[stencil.index[i]];
				if (stencil.index[i] == index)
					weight += stencil.coefficient[i];
			}
			d_energy += prefactor * (-cos(kPi * (old_sum + weight * delta)) + cos(kPi * old_sum));
		};

	const double e_pref = 1.0 / KU;
	const double b_pref = 1.0;

	if (index >= nSites) // A^0(s): appears in E^xy at s+-x+-y and in E^xx at s+-x, s+-y
	{
		// Given a A0 index, find the corresponding type1 index it sits over
		// Since There are twice as many indices (type 1 and 2) than A0 indices
		// One has to first multiply by two and afterwards offset towards the type 1
		// site, the offset depending on the current row
		const int pair = 2 * (index - nSites);
		const int s = pair + lattice.coord(pair, 1) % 2;

		for (int site : { shifted(s, +1, +1, 0), shifted(s, -1, +1, 0),
						  shifted(s, -1, -1, 0), shifted(s, +1, -1, 0) })
			accumulate(electricStencil_xy(site), e_pref);

		for (int site : { shifted(s, +1, 0, 0), shifted(s, -1, 0, 0),
						  shifted(s, 0, +1, 0), shifted(s, 0, -1, 0) })
			accumulate(electricStencil_xx(site), e_pref);
	}
	else if (isType1(index)) // A^xy(s): appears in E^xy at s, s-t and in B at s+-x, s+-y
	{
		const int s = index;

		for (int site : { s, shifted(s, 0, 0, -1) })
			accumulate(electricStencil_xy(site), e_pref);

		for (int site : { shifted(s, +1, 0, 0), shifted(s, -1, 0, 0),
						  shifted(s, 0, +1, 0), shifted(s, 0, -1, 0) })
			accumulate(magneticStencil(site), b_pref);
	}
	else // A^xx(s): appears in E^xx at s, s-t and in B at s+-x+-y
	{
		const int s = index;

		for (int site : { s, shifted(s, 0, 0, -1) })
			accumulate(electricStencil_xx(site), e_pref);

		for (int site : { shifted(s, +1, +1, 0), shifted(s, -1, +1, 0),
						  shifted(s, -1, -1, 0), shifted(s, +1, -1, 0) })
			accumulate(magneticStencil(site), b_pref);
	}

	return d_energy;
}

// A gauge transformation with f a delta function: f = random_shift on one type 1 site s, 0 elsewhere.
//   A^0(r)  --> A^0(r)  + f(r+t) - f(r)
//   A^xy(r) --> A^xy(r) + f(r+x+y) - f(r-x+y) + f(r-x-y) - f(r+x-y)
//   A^xx(r) --> A^xx(r) + f(r+x) - f(r+y) + f(r-x) - f(r-y)
// Tests show that gauge transformations are useless as overrelaxations, i.e. do not help with exploring
// The hilbert space and do not help with thermalization.
void Spiderweb_corrected::gaugeTransformation(int index, std::mt19937& rng)
{
	std::uniform_real_distribution<double> shift_dst(-1.0, 1.0);

	// f lives on type 1 sites. Every pair (2k, 2k+1) holds one, and three variables belong to it
	// (its A^0, its A^xy and the A^xx of its type 2 partner), so a uniformly drawn variable gives a
	// uniformly drawn type 1 site.
	const int pair = 2 * (index >= nSites ? index - nSites : index / 2);
	const int s = pair + lattice.coord(pair, 1) % 2;

	double random_shift = shift_dst(rng);

	// A^0(s) --> A^0(s) + f(s+t) - f(s) = A^0(s) - f(s)
	fields[a0Index(s)] += -random_shift;
	// A^0(s-t) --> A^0(s-t) + f(s) - f(s-t) = A^0(s-t) + f(s)
	fields[a0Index(shifted(s, 0, 0, -1))] += random_shift;

	// A^xy on the four diagonal type 1 neighbours; f(s) enters A^xy(r) as the term f(r + (s - r)).
	fields[tensorIndex(shifted(s, -1, -1, 0))] += random_shift;  // + f(r+x+y)
	fields[tensorIndex(shifted(s, +1, -1, 0))] += -random_shift; // - f(r-x+y)
	fields[tensorIndex(shifted(s, +1, +1, 0))] += random_shift;  // + f(r-x-y)
	fields[tensorIndex(shifted(s, -1, +1, 0))] += -random_shift; // - f(r+x-y)

	// A^xx on the four nearest type 2 neighbours, in the same way.
	fields[tensorIndex(shifted(s, -1, 0, 0))] += random_shift;   // + f(r+x)
	fields[tensorIndex(shifted(s, 0, -1, 0))] += -random_shift;  // - f(r+y)
	fields[tensorIndex(shifted(s, +1, 0, 0))] += random_shift;   // + f(r-x)
	fields[tensorIndex(shifted(s, 0, +1, 0))] += -random_shift;  // - f(r-y)
}

std::vector<std::string> Spiderweb_corrected::observableNames() const
{
	return { "Energy", "HelicityB", "HelicityExx", "HelicityExy", "InstantonDensity" };
}

double Spiderweb_corrected::instantonDensity() const
{
	int n_instantons = 0;
	for (int m : getInstantons())
		n_instantons += std::abs(m);

	// Instantons live on type 2 sites only, so there are nSites / 2 places for one.
	return n_instantons / (nSites / 2.0);
}

// Dont trust the electric helicity modulus (yet). One could probably also implement a genuine
// helicity modulus as a twist in the A fields along temporal directions. Still have to think about that
std::array<double, 3> Spiderweb_corrected::helicitySamples(double temperature) const
{
	// Sum of cos and of sin over every term of one kind: B, E^xx, E^xy.
	double cos_b = 0.0, sin_b = 0.0, cos_xx = 0.0, sin_xx = 0.0, cos_xy = 0.0, sin_xy = 0.0;
	for (int site = 0; site < nSites; site++)
	{
		if (isType1(site))
		{
			const double e = kPi * electricField_xy(site);
			cos_xy += cos(e); sin_xy += sin(e);
		}
		else
		{
			const double b = kPi * magneticField(site);
			const double e = kPi * electricField_xx(site);
			cos_b += cos(b); sin_b += sin(b);
			cos_xx += cos(e); sin_xx += sin(e);
		}
	}

	// Each kind of term sits on one sublattice, so there are nSites / 2 of each.
	const double n_terms = nSites / 2.0;
	return { (cos_b - sin_b * sin_b / temperature) / n_terms,
			 (cos_xx - sin_xx * sin_xx / (temperature * KU)) / n_terms,
			 (cos_xy - sin_xy * sin_xy / (temperature * KU)) / n_terms };
}

std::vector<double> Spiderweb_corrected::measure(double temperature) const
{
	const std::array<double, 3> helicity = helicitySamples(temperature);
	return { getEnergy(), helicity[0], helicity[1], helicity[2], instantonDensity() };
}

std::vector<Channel> Spiderweb_corrected::channels() const
{
	// Every range is known from the model alone, so the colours never depend on the state shown.
	// The fluxes are wrapped onto [-1, 1), in units of pi. The local energy is drawn from each
	// site's own minimum, so the ground state is uniformly black on both sublattices; the largest
	// value is then that of a type 2 site with B and E^xx both at pi. E^xy and E^xx live on
	// different sublattices, so one channel shows both: E^xy on type 1 sites, E^xx on type 2.
	// The same goes for A: A^xy on type 1 sites, A^xx on type 2.
	return { { "local energy", ChannelKind::Magnitude, 0.0, 2.0 + 2.0 / KU },
			 { "E", ChannelKind::Signed, -1.0, 1.0 },
			 { "B", ChannelKind::Signed, -1.0, 1.0 },
			 { "A_0", ChannelKind::Angle },
			 { "A", ChannelKind::Angle },
			 { "instantons", ChannelKind::Integer } };
}

void Spiderweb_corrected::fillChannel(int channel, std::vector<double>& out) const
{
	// A field is only defined on its own sublattice; the other one is left without a value.
	const double no_value = std::numeric_limits<double>::quiet_NaN();
	out.assign(lattice.size(), no_value);

	if (channel == 0)
	{
		out = getLocalEnergies();
		for (int site = 0; site < nSites; site++)
			out[site] += isType1(site) ? 1.0 / KU : 1.0 + 1.0 / KU;
		return;
	}
	if (channel == 5)
	{
		// Violations of the Bianchi identity, see getInstantons(). They live on type 2 sites.
		const std::vector<int> instantons = getInstantons();
		for (int site = 0; site < nSites; site++)
			if (!isType1(site))
				out[site] = (double)instantons[site];
		return;
	}

	for (int site = 0; site < nSites; site++)
	{
		const bool type1 = isType1(site);
		switch (channel)
		{
		case 1: out[site] = mcf::mapToCircle(type1 ? electricField_xy(site) : electricField_xx(site)); break;
		case 2: if (!type1) out[site] = mcf::mapToCircle(magneticField(site)); break;
		case 3: if (type1)  out[site] = fields[a0Index(site)]; break;
		case 4: out[site] = fields[tensorIndex(site)]; break;
		default: break;
		}
	}
}
