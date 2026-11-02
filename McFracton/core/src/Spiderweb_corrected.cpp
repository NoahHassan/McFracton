#include "Spiderweb_corrected.h"

#include <assert.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
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
	buildPatterns();
}

void Spiderweb_corrected::buildPatterns()
{
	// Sublattice and neighbours first: isType1() and shifted() read them, and the stencil functions
	// used below call both.
	type1.resize(nSites);
	for (int site = 0; site < nSites; site++)
		type1[site] = (lattice.coord(site, 0) + lattice.coord(site, 1)) % 2 == 0;
	neighbors.resize(6 * (size_t)nSites);
	for (int site = 0; site < nSites; site++)
		for (int axis = 0; axis < 3; axis++)
		{
			neighbors[6 * site + 2 * axis + 0] = lattice.neighbor(site, axis, -1);
			neighbors[6 * site + 2 * axis + 1] = lattice.neighbor(site, axis, +1);
		}

	// Offset of a variable from a site, as the shortest periodic displacement. On the smallest
	// lattices +2 and -2 (or +1 and -1 in t) are the same site, and then either is right.
	auto displacement = [&](int from, int to, int axis) {
		const int L = lattice.extent(axis);
		int d = ((lattice.coord(to, axis) - lattice.coord(from, axis)) % L + L) % L;
		if (2 * d >= L)
			d -= L;
		return d;
	};
	auto makeTerm = [&](int kind, int term_site, int around, int variable) {
		const Stencil s = stencil(kind, term_site);
		Term term;
		term.kind = kind;
		term.size = s.size;
		for (int i = 0; i < s.size; i++)
		{
			const int site = variableSite(s.index[i]);
			Tap& tap = term.taps[i];
			tap.dx = displacement(around, site, 0);
			tap.dy = displacement(around, site, 1);
			tap.dt = displacement(around, site, 2);
			tap.a0 = s.index[i] >= nSites ? 1 : 0;
			tap.coefficient = s.coefficient[i];
			// The weight exactly as proposeUpdate has always worked it out: the sum of the
			// coefficients the updated variable has in the term.
			if (s.index[i] == variable)
				term.weight += s.coefficient[i];
		}
		return term;
	};

	// The patterns, read off one variable of each kind (A^0 of site 0, A^xy of site 0, A^xx of
	// site 1) and one term of each kind (B and E^xx at site 1, E^xy at site 0).
	const int reference_variable[3] = { nSites, 0, 1 };
	for (int c = 0; c < 3; c++)
	{
		const int index = reference_variable[c];
		std::array<TermSite, 8> terms;
		Pattern& pattern = variable_patterns[c];
		pattern.size = couplingTerms(index, terms);
		for (int k = 0; k < pattern.size; k++)
			pattern.terms[k] = makeTerm(terms[k].kind, terms[k].site, variableSite(index), index);
	}
	term_patterns[kMagnetic] = makeTerm(kMagnetic, 1, 1, -1);
	term_patterns[kElectricXX] = makeTerm(kElectricXX, 1, 1, -1);
	term_patterns[kElectricXY] = makeTerm(kElectricXY, 0, 0, -1);

	// Wrap tables, as wide as the patterns reach.
	reach_xy = 0;
	reach_t = 0;
	auto widen = [&](const Term& term) {
		for (int i = 0; i < term.size; i++)
		{
			reach_xy = std::max({ reach_xy, std::abs(term.taps[i].dx), std::abs(term.taps[i].dy) });
			reach_t = std::max(reach_t, std::abs(term.taps[i].dt));
		}
	};
	for (const Pattern& pattern : variable_patterns)
		for (int k = 0; k < pattern.size; k++)
			widen(pattern.terms[k]);
	for (const Term& term : term_patterns)
		widen(term);
	const int L = linear_size, T = temporal_size;
	wrap_x.resize(L + 2 * reach_xy);
	wrap_y.resize(L + 2 * reach_xy);
	wrap_t.resize(T + 2 * reach_t);
	for (int i = 0; i < (int)wrap_x.size(); i++)
	{
		wrap_x[i] = ((i - reach_xy) % L + L) % L;
		wrap_y[i] = wrap_x[i] * L;
	}
	for (int i = 0; i < (int)wrap_t.size(); i++)
		wrap_t[i] = ((i - reach_t) % T + T) % T * L * L;

	// Division by L as a multiplication by 2^40 / L, rounded up. It is exact as long as n * L < 2^40,
	// far beyond any lattice that fits in memory; the loop below checks it for every n it is used on.
	l_reciprocal = ((std::uint64_t)1 << 40) / (std::uint64_t)L + 1;

	// Every pattern is checked against the stencil functions, for every variable and every site:
	// same variables, same coefficients, same order, same weights. Anything that is not exactly
	// a shifted copy of the reference refuses to build.
	auto check = [](bool ok, const char* what) {
		if (!ok)
			throw std::logic_error(std::string("Spiderweb_corrected: stencil pattern check failed: ") + what);
	};
	for (int n = 0; n < nSites; n++)
		check(divideByL(n) == n / L, "division by L");
	auto tapIndex = [&](const Tap& tap, const Origin& o) {
		const int site = o.x[tap.dx] + o.y[tap.dy] + o.t[tap.dt];
		return tap.a0 ? a0Index(site) : tensorIndex(site);
	};
	auto sameTerm = [&](const Term& term, const Stencil& s, const Origin& o, int variable) {
		// tapSum unrolls each kind of term to its fixed size.
		if (term.size != s.size || term.size != (term.kind == kMagnetic ? kMagneticSize : kElectricSize))
			return false;
		double weight = 0.0;
		for (int i = 0; i < s.size; i++)
		{
			if (tapIndex(term.taps[i], o) != s.index[i] || term.taps[i].coefficient != s.coefficient[i])
				return false;
			// proposeUpdate relies on every variable entering a term at most once.
			for (int j = 0; j < i; j++)
				if (s.index[j] == s.index[i])
					return false;
			if (s.index[i] == variable)
				weight += s.coefficient[i];
		}
		return variable < 0 || (weight == term.weight && (weight == 1.0 || weight == -1.0));
	};
	std::vector<int> appearances(numVariables(), 0);
	for (int site = 0; site < nSites; site++)
	{
		const Origin o = origin(site);
		check(o.x == wrap_x.data() + reach_xy + lattice.coord(site, 0)
			&& o.y == wrap_y.data() + reach_xy + lattice.coord(site, 1)
			&& o.t == wrap_t.data() + reach_t + lattice.coord(site, 2), "site coordinates");
		for (int kind : { kMagnetic, kElectricXX, kElectricXY })
		{
			if ((kind == kElectricXY) != isType1(site))
				continue;
			const Stencil s = stencil(kind, site);
			check(sameTerm(term_patterns[kind], s, o, -1), "term");
			for (int i = 0; i < s.size; i++)
				appearances[s.index[i]]++;
		}
	}
	for (int index = 0; index < numVariables(); index++)
	{
		// The same lookup proposeUpdate makes, checked against the slow, explicit one.
		int kind;
		const Origin o = variableOrigin(index, kind);
		const Origin expected = origin(variableSite(index));
		check(kind == (index >= nSites ? 0 : isType1(index) ? 1 : 2)
			&& o.x == expected.x && o.y == expected.y && o.t == expected.t, "variable site");
		const Pattern& pattern = variable_patterns[kind];
		std::array<TermSite, 8> terms;
		const int n_terms = couplingTerms(index, terms);
		check(n_terms == pattern.size, "number of terms");
		for (int k = 0; k < n_terms; k++)
			check(pattern.terms[k].kind == terms[k].kind
				&& sameTerm(pattern.terms[k], stencil(terms[k].kind, terms[k].site), o, index), "coupling");
		// And nothing is missed: a variable enters exactly the terms it is coupled to.
		check(appearances[index] == n_terms, "a variable enters a term it is not coupled to");
	}
}

Spiderweb_corrected::Stencil Spiderweb_corrected::stencil(int kind, int site_index) const
{
	switch (kind)
	{
	case kMagnetic: return magneticStencil(site_index);
	case kElectricXX: return electricStencil_xx(site_index);
	default: return electricStencil_xy(site_index);
	}
}

int Spiderweb_corrected::variableSite(int index) const
{
	if (index < nSites)
		return index;
	// Every pair of sites (2k, 2k+1) holds one type 1 site, the even one in an even row and the
	// odd one in an odd row (x runs fastest and linear_size is even).
	const int pair = 2 * (index - nSites);
	return pair + lattice.coord(pair, 1) % 2;
}

int Spiderweb_corrected::couplingTerms(int index, std::array<TermSite, 8>& terms) const
{
	const int s = variableSite(index);
	int n = 0;
	auto add = [&](int kind, int site) { terms[n++] = { kind, site }; };

	if (index >= nSites) // A^0(s): appears in E^xy at s+-x+-y and in E^xx at s+-x, s+-y
	{
		for (int site : { shifted(s, +1, +1, 0), shifted(s, -1, +1, 0),
						  shifted(s, -1, -1, 0), shifted(s, +1, -1, 0) })
			add(kElectricXY, site);
		for (int site : { shifted(s, +1, 0, 0), shifted(s, -1, 0, 0),
						  shifted(s, 0, +1, 0), shifted(s, 0, -1, 0) })
			add(kElectricXX, site);
	}
	else if (isType1(index)) // A^xy(s): appears in E^xy at s, s-t and in B at s+-x, s+-y
	{
		for (int site : { s, shifted(s, 0, 0, -1) })
			add(kElectricXY, site);
		for (int site : { shifted(s, +1, 0, 0), shifted(s, -1, 0, 0),
						  shifted(s, 0, +1, 0), shifted(s, 0, -1, 0) })
			add(kMagnetic, site);
	}
	else // A^xx(s): appears in E^xx at s, s-t and in B at s+-x+-y
	{
		for (int site : { s, shifted(s, 0, 0, -1) })
			add(kElectricXX, site);
		for (int site : { shifted(s, +1, +1, 0), shifted(s, -1, +1, 0),
						  shifted(s, -1, -1, 0), shifted(s, +1, -1, 0) })
			add(kMagnetic, site);
	}
	return n;
}

Spiderweb_corrected::Origin Spiderweb_corrected::origin(int site_index) const
{
	const int row = divideByL(site_index);
	const int t = divideByL(row);
	const int x = site_index - row * linear_size;
	const int y = row - t * linear_size;
	return { wrap_x.data() + reach_xy + x, wrap_y.data() + reach_xy + y, wrap_t.data() + reach_t + t };
}

// variableSite() and the pattern a variable follows (0: A^0, 1: A^xy, 2: A^xx) in one go, from a
// single set of coordinates; proposeUpdate calls it once per update.
Spiderweb_corrected::Origin Spiderweb_corrected::variableOrigin(int index, int& pattern) const
{
	const bool is_a0 = index >= nSites;
	// A^0 sits on the type 1 site of the pair (2k, 2k+1): the even one in an even row, the odd one
	// in an odd row.
	const int first = is_a0 ? 2 * (index - nSites) : index;
	const int row = divideByL(first);
	const int t = divideByL(row);
	const int y = row - t * linear_size;
	int x = first - row * linear_size;
	if (is_a0)
		x += y % 2;
	pattern = is_a0 ? 0 : (x + y) % 2 == 0 ? 1 : 2;
	return { wrap_x.data() + reach_xy + x, wrap_y.data() + reach_xy + y, wrap_t.data() + reach_t + t };
}

// Sum of coefficient * field over one stencil, in stencil order. It is the sum the stencils have
// always been evaluated with: same variables, same order, same bits. The size is a template
// argument only so that the compiler can unroll the loop; buildPatterns checks that every B term
// has kMagneticSize taps and every E term kElectricSize.
template <int N>
double Spiderweb_corrected::tapSum(const Term& term, const Origin& o) const
{
	double total = 0.0;
	for (int i = 0; i < N; i++)
	{
		const Tap& tap = term.taps[i];
		const int site = o.x[tap.dx] + o.y[tap.dy] + o.t[tap.dt];
		total += tap.coefficient * fields[tap.a0 ? nSites + site / 2 : site];
	}
	return total;
}

double Spiderweb_corrected::termSum(const Term& term, const Origin& o) const
{
	return term.kind == kMagnetic ? tapSum<kMagneticSize>(term, o) : tapSum<kElectricSize>(term, o);
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

// Type 1 when x + y is even, filled in by buildPatterns().
bool Spiderweb_corrected::isType1(int site_index) const
{
	return type1[site_index] != 0;
}

// To avoid having two fields on type 1 (A^xy and A^0) and one field on type 2 (A^xx)
// A layout is chosen by which first all spatial components are listed and the
// A0 fields come after in the array
int Spiderweb_corrected::a0Index(int site_index) const
{
	assert(isType1(site_index)); // Important because A^0 lives on links connecting crossed sties
	return nSites + site_index / 2;
}

// Every step is at most one site along each axis, so the neighbour table covers all of them.
int Spiderweb_corrected::shifted(int site_index, int dx, int dy, int dt) const
{
	assert(-1 <= dx && dx <= 1 && -1 <= dy && dy <= 1 && -1 <= dt && dt <= 1);
	int site = site_index;
	if (dx != 0) site = neighbors[6 * site + 0 + (dx > 0)];
	if (dy != 0) site = neighbors[6 * site + 2 + (dy > 0)];
	if (dt != 0) site = neighbors[6 * site + 4 + (dt > 0)];
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

double Spiderweb_corrected::magneticField(int site_index) const
{
	assert(!isType1(site_index));
	return termSum(term_patterns[kMagnetic], origin(site_index));
}

double Spiderweb_corrected::electricField_xx(int site_index) const
{
	assert(!isType1(site_index));
	return termSum(term_patterns[kElectricXX], origin(site_index));
}

double Spiderweb_corrected::electricField_xy(int site_index) const
{
	assert(isType1(site_index));
	return termSum(term_patterns[kElectricXY], origin(site_index));
}

// The Hamiltonian per site, from the cosines of the terms on it:
//   type 1: - 1/KU cos(E^xy),   type 2: - cos(B) - 1/KU cos(E^xx)
double Spiderweb_corrected::type1Energy(double cos_xy) const
{
	return -cos_xy / KU;
}

double Spiderweb_corrected::type2Energy(double cos_b, double cos_xx) const
{
	return -cos_b - cos_xx / KU;
}

double Spiderweb_corrected::localEnergy(int site_index) const
{
	if (isType1(site_index))
		return type1Energy(cos(kPi * electricField_xy(site_index)));
	return type2Energy(cos(kPi * magneticField(site_index)), cos(kPi * electricField_xx(site_index)));
}

std::vector<double> Spiderweb_corrected::getLocalEnergies() const
{
	std::vector<double> localEnergies(nSites);
	for (int site = 0; site < nSites; site++)
		localEnergies[site] = localEnergy(site);
	return localEnergies;
}

// One pass over the lattice that computes every field strength and its cosine once, and from them
// whatever the observables need. Each sum is taken site by site in index order, which is the order
// the observables have always been summed in, so measure() gives the same bits as the separate
// getEnergy(), helicitySamples() and instantonDensity() did.
void Spiderweb_corrected::latticePass(LatticeSums& sums, bool helicity, bool wrapped) const
{
	sums = LatticeSums{};
	if (wrapped)
	{
		// B and E wrapped onto the circle, for the instantons.
		sums.b.assign(nSites, 0.0);
		sums.e.assign(nSites, 0.0);
	}
	for (int site = 0; site < nSites; site++)
	{
		if (isType1(site))
		{
			const double exy = electricField_xy(site);
			const double cos_xy = cos(kPi * exy);
			sums.energy += type1Energy(cos_xy);
			if (helicity)
			{
				sums.cos_xy += cos_xy;
				sums.sin_xy += sin(kPi * exy);
			}
			if (wrapped)
				sums.e[site] = mcf::mapToCircle(exy);
		}
		else
		{
			const double b = magneticField(site);
			const double exx = electricField_xx(site);
			const double cos_b = cos(kPi * b);
			const double cos_xx = cos(kPi * exx);
			sums.energy += type2Energy(cos_b, cos_xx);
			if (helicity)
			{
				sums.cos_b += cos_b;
				sums.sin_b += sin(kPi * b);
				sums.cos_xx += cos_xx;
				sums.sin_xx += sin(kPi * exx);
			}
			if (wrapped)
			{
				sums.b[site] = mcf::mapToCircle(b);
				sums.e[site] = mcf::mapToCircle(exx);
			}
		}
	}
}

// The sum of getLocalEnergies(), in the same order, without building the vector: McMachine calls
// this once per sweep while measuring.
double Spiderweb_corrected::getEnergy() const
{
	LatticeSums sums;
	latticePass(sums, false, false);
	return sums.energy;
}

// Instantons are computed as violations of the gauss law Dt B + Dxx E^xy - Dxy E^xx = 0, from B and
// E wrapped onto the circle. r has to be a type 2 site.
int Spiderweb_corrected::instantonAt(int r, const LatticeSums& sums) const
{
	const std::vector<double>& b = sums.b;
	const std::vector<double>& e = sums.e;
	const double wrapped_gauss = b[shifted(r, 0, 0, +1)] - b[r]
		// Dxx [E^xy]: the four type 1 neighbours
		+ e[shifted(r, +1, 0, 0)] - e[shifted(r, 0, +1, 0)] + e[shifted(r, -1, 0, 0)] - e[shifted(r, 0, -1, 0)]
		// - Dxy [E^xx]: the four diagonal type 2 neighbours
		- (e[shifted(r, +1, +1, 0)] - e[shifted(r, -1, +1, 0)] + e[shifted(r, -1, -1, 0)] - e[shifted(r, +1, -1, 0)]);
	return (int)std::lround(wrapped_gauss / 2.0);
}

std::vector<int> Spiderweb_corrected::getInstantons() const
{
	LatticeSums sums;
	latticePass(sums, false, true);

	std::vector<int> instantons(nSites, 0);
	for (int r = 0; r < nSites; r++)
		if (!isType1(r))
			instantons[r] = instantonAt(r, sums);
	return instantons;
}

double Spiderweb_corrected::proposeUpdate(int index, double delta) const
{
	double d_energy = 0.0;

	// The energy difference of one cosine term, given its current sum and the coefficient (weight)
	// of the variable at index in it.
	auto accumulate = [&](double old_sum, double weight, double prefactor)
		{
			d_energy += prefactor * (-cos(kPi * (old_sum + weight * delta)) + cos(kPi * old_sum));
		};

	const double e_pref = 1.0 / KU;
	const double b_pref = 1.0;

	// The terms the variable enters, in the order of couplingTerms(), which is the order this
	// function has always added them up in:
	//   A^0(s):  E^xy at s+-x+-y, then E^xx at s+-x, s+-y
	//   A^xy(s): E^xy at s, s-t,  then B at s+-x, s+-y
	//   A^xx(s): E^xx at s, s-t,  then B at s+-x+-y
	int kind;
	const Origin o = variableOrigin(index, kind);
	const Pattern& pattern = variable_patterns[kind];
	for (int k = 0; k < pattern.size; k++)
	{
		const Term& term = pattern.terms[k];
		accumulate(termSum(term, o), term.weight, term.kind == kMagnetic ? b_pref : e_pref);
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

double Spiderweb_corrected::instantonDensity(const LatticeSums& sums) const
{
	int n_instantons = 0;
	for (int r = 0; r < nSites; r++)
		if (!isType1(r))
			n_instantons += std::abs(instantonAt(r, sums));

	// Instantons live on type 2 sites only, so there are nSites / 2 places for one.
	return n_instantons / (nSites / 2.0);
}

double Spiderweb_corrected::instantonDensity() const
{
	LatticeSums sums;
	latticePass(sums, false, true);
	return instantonDensity(sums);
}

// Dont trust the electric helicity modulus (yet). One could probably also implement a genuine
// helicity modulus as a twist in the A fields along temporal directions. Still have to think about that
std::array<double, 3> Spiderweb_corrected::helicitySamples(const LatticeSums& sums, double temperature) const
{
	// sums holds the sum of cos and of sin over every term of one kind: B, E^xx, E^xy.
	// Each kind of term sits on one sublattice, so there are nSites / 2 of each.
	const double n_terms = nSites / 2.0;
	return { (sums.cos_b - sums.sin_b * sums.sin_b / temperature) / n_terms,
			 (sums.cos_xx - sums.sin_xx * sums.sin_xx / (temperature * KU)) / n_terms,
			 (sums.cos_xy - sums.sin_xy * sums.sin_xy / (temperature * KU)) / n_terms };
}

std::array<double, 3> Spiderweb_corrected::helicitySamples(double temperature) const
{
	LatticeSums sums;
	latticePass(sums, true, false);
	return helicitySamples(sums, temperature);
}

// All observables from a single pass over the lattice; see latticePass.
std::vector<double> Spiderweb_corrected::measure(double temperature) const
{
	LatticeSums sums;
	latticePass(sums, true, true);
	const std::array<double, 3> helicity = helicitySamples(sums, temperature);
	return { sums.energy, helicity[0], helicity[1], helicity[2], instantonDensity(sums) };
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
