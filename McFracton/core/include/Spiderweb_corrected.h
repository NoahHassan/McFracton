#pragma once

#include <array>
#include <random>

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
	double sum(const Stencil& stencil) const;
	int shifted(int site_index, int dx, int dy, int dt) const;
};
