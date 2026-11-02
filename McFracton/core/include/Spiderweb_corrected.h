#pragma once

#include <array>
#include <cstdint>
#include <random>
#include <vector>

#include "Lattice.h"
#include "System.h"

// The spiderweb model on a square lattice with a two-site basis. Axes: x, y, t.
//
//   H = - sum_type2 ( cos(B) + 1/KU cos(E^xx) ) - sum_type1 1/KU cos(E^xy)
//
// The basis is the checkerboard of the xy plane: a site (x, y, t) is type 1 when x + y is even and
// type 2 when it is odd, for every t. That is the assignment under which every field in the
// definitions below sits on the sublattice it is said to live on, and it needs an even linear size
// to close periodically.
//
//   type 1 carries A^0 and A^xy, and the term E^xy
//   type 2 carries A^xx,          and the terms B and E^xx
//
//   B(r)    = A^xy(r+x) - A^xx(r+x+y) - A^xy(r+y) + A^xx(r-x+y)
//           + A^xy(r-x) - A^xx(r-x-y) - A^xy(r-y) + A^xx(r+x-y)
//   E^xy(r) = A^0(r+x+y) - A^0(r-x+y) + A^0(r-x-y) - A^0(r+x-y) - A^xy(r+t) + A^xy(r)
//   E^xx(r) = A^0(r+x)   - A^0(r+y)   + A^0(r-x)   - A^0(r-y)   - A^xx(r+t) + A^xx(r)
//
// Fields are in units of pi, like everywhere else in the code, so cos(B) is cos(pi * B) in terms
// of the stored numbers and every field has period 2.
//
// Variables: there are 3/2 per site. The first nSites are the tensor component living on that
// site (A^xy on type 1, A^xx on type 2), so their index is the site index. The remaining nSites/2
// are A^0, one per type 1 site.
//
// Gauge transformations, under which all three are invariant, with f on the type 1 sites:
//   A^0(r)  --> A^0(r)  + f(r+t) - f(r)
//   A^xy(r) --> A^xy(r) + f(r+x+y) - f(r-x+y) + f(r-x-y) - f(r+x-y)
//   A^xx(r) --> A^xx(r) + f(r+x) - f(r+y) + f(r-x) - f(r-y)
// gaugeTransformation applies one with f nonzero on a single site. It is not an update of the
// Monte Carlo: overrelaxation is deliberately not implemented, so it falls back to the base class throw.
class Spiderweb_corrected : public System {
public:
	Spiderweb_corrected(int linear_size, int temporal_size, double KU);
	~Spiderweb_corrected() override = default;
public:
	void randomize(std::mt19937& rng) override;
	void setGroundState() override;
	double getEnergy() const override;
	std::vector<double> getLocalEnergies() const;
	double localEnergy(int site_index) const;
	double proposeUpdate(int index, double delta) const override;
	std::vector<std::string> observableNames() const override;
	std::vector<double> measure(double temperature) const override;
	std::vector<Channel> channels() const override;
	void fillChannel(int channel, std::vector<double>& out) const override;
public:
	// Gauge transformation with f random on the type 1 site belonging to variable `index`.
	void gaugeTransformation(int index, std::mt19937& rng);
	bool isType1(int site_index) const;
	// Variable index of the tensor component on a site: A^xy on type 1, A^xx on type 2.
	int tensorIndex(int site_index) const { return site_index; }
	// Variable index of A^0 on a type 1 site.
	int a0Index(int site_index) const;
	double magneticField(int site_index) const;
	double electricField_xx(int site_index) const;
	double electricField_xy(int site_index) const;
	std::vector<int> getInstantons() const;
	double instantonDensity() const;
	std::array<double, 3> helicitySamples(double temperature) const;
	// The fields entering one cosine, as (variable index, coefficient) pairs.
	struct Stencil {
		std::array<int, 8> index{};
		std::array<double, 8> coefficient{};
		int size = 0;
		void add(int variable, double factor) { index[size] = variable; coefficient[size] = factor; size++; }
	};
	Stencil magneticStencil(int site_index) const;
	Stencil electricStencil_xx(int site_index) const;
	Stencil electricStencil_xy(int site_index) const;
public:
	const int linear_size;
	const int temporal_size;
	const int nSites;
	const double KU;
private:
	enum TermKind { kMagnetic = 0, kElectricXX = 1, kElectricXY = 2 };
	struct TermSite { int kind; int site; };
	int shifted(int site_index, int dx, int dy, int dt) const;
	Stencil stencil(int kind, int site_index) const;
	// The site a variable sits on: its own for A^xy and A^xx, the type 1 site of its pair for A^0.
	int variableSite(int index) const;
	// The terms variable `index` enters, in the order proposeUpdate adds them up:
	//   A^0:  4 x E^xy (s+x+y, s-x+y, s-x-y, s+x-y), then 4 x E^xx (s+x, s-x, s+y, s-y)
	//   A^xy: 2 x E^xy (s, s-t),                      then 4 x B    (s+x, s-x, s+y, s-y)
	//   A^xx: 2 x E^xx (s, s-t),                      then 4 x B    (s+x+y, s-x+y, s-x-y, s+x-y)
	int couplingTerms(int index, std::array<TermSite, 8>& terms) const;
	void buildPatterns();

	// What one pass over the lattice gives the observables: the energy, the sums of cos and sin of
	// every kind of term (helicity moduli) and, if asked for, B and E wrapped onto the circle
	// (instantons).
	struct LatticeSums {
		double energy = 0.0;
		double cos_b = 0.0, sin_b = 0.0, cos_xx = 0.0, sin_xx = 0.0, cos_xy = 0.0, sin_xy = 0.0;
		std::vector<double> b, e;
	};
	void latticePass(LatticeSums& sums, bool helicity, bool wrapped) const;
	double type1Energy(double cos_xy) const;
	double type2Energy(double cos_b, double cos_xx) const;
	int instantonAt(int type2_site, const LatticeSums& sums) const;
	double instantonDensity(const LatticeSums& sums) const;
	std::array<double, 3> helicitySamples(const LatticeSums& sums, double temperature) const;
private:
	// Stencil patterns, built once in the constructor from the stencil functions above, so the
	// model is still defined in exactly one place, and checked against them for every variable and
	// every site. They only replace index arithmetic: every sum runs over the same variables, with
	// the same coefficients, in the same order as before, so all results are bit-identical to
	// building the stencils on the fly.
	//
	// The model is translation invariant, so the terms of one kind of variable (or the stencil of
	// one kind of term) are the same everywhere up to a shift. A pattern stores them as offsets
	// (dx, dy, dt) from the variable's site, and the index of a field is found from the site's
	// coordinates through the small wrap tables below. Unlike a table per variable this needs no
	// memory that grows with the lattice, so an update reads nothing but the fields themselves.
	struct Tap {
		int dx = 0, dy = 0, dt = 0;
		int a0 = 0;          // 1: the A^0 of that site, 0: its tensor component
		double coefficient = 0.0;
	};
	struct Term {
		int kind = kMagnetic;
		int size = 0;
		double weight = 0.0; // coefficient of the updated variable in this term
		std::array<Tap, 8> taps{};
	};
	struct Pattern {
		int size = 0;
		std::array<Term, 8> terms{};
	};
	// Where wrap tables have to be read for a site: x, y and t already folded in.
	struct Origin { const int* x; const int* y; const int* t; };
	Origin origin(int site_index) const;
	Origin variableOrigin(int index, int& pattern) const;
	static constexpr int kMagneticSize = 8;
	static constexpr int kElectricSize = 6;
	template <int N>
	double tapSum(const Term& term, const Origin& o) const;
	double termSum(const Term& term, const Origin& o) const;
	// n / linear_size as a multiplication, for 0 <= n < nSites; buildPatterns checks it for all n.
	int divideByL(int n) const { return int((std::uint64_t(n) * l_reciprocal) >> 40); }
	std::uint64_t l_reciprocal = 0;

	std::array<Pattern, 3> variable_patterns; // A^0, A^xy, A^xx: the terms each one enters
	std::array<Term, 3> term_patterns;        // B, E^xx, E^xy around their own site
	// wrap_x[reach_xy + x + dx] is (x + dx) mod L, wrap_y likewise times L, wrap_t times L^2.
	int reach_xy = 0, reach_t = 0;
	std::vector<int> wrap_x, wrap_y, wrap_t;
	// Neighbour of every site one step along each axis: neighbors[6 * site + 2 * axis + (step > 0)].
	std::vector<int> neighbors;
	std::vector<unsigned char> type1;
};
