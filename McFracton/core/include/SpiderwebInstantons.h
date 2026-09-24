#pragma once

// Instanton diagnostics for Spiderweb, written against its public interface only. They back the
// "instantons m" and "Z4 steps" channels and change nothing about the physics.
//
// Fluxes per site r (code units: angles in units of pi, period 2):
//   e_xx(r) = Q_xx A0 - D0 A_xx,   e_xy(r) = Q_xy A0 - D0 A_xy,   b(r) = Q_xy A_xx - Q_xx A_xy
//   with Q_xx = 4 Dx Dy, Q_xy = Dx^2 - Dy^2 and forward differences, exactly as in Spiderweb.cpp.
//
// Bianchi identity (exact for the unwrapped sums):  D0 b + Q_xy e_xx - Q_xx e_xy = 0.
// With wrapped fluxes [f] = mcf::mapToCircle(f) it equals 2 m(r) with m integer: the instanton
// density, the rank-2 counterpart of the monopole divergence in AbelianGaugeSquare::getMonopoles.
//
//   2 m(x,y,t) = [b](x,y,t+1) - [b](x,y,t)
//              + [e_xx](x+2,y,t) - 2[e_xx](x+1,y,t) - [e_xx](x,y+2,t) + 2[e_xx](x,y+1,t)
//              - 4 ( [e_xy](x+1,y+1,t) - [e_xy](x+1,y,t) - [e_xy](x,y+1,t) + [e_xy](x,y,t) )
//
// Z4 steps: A_xy enters b only as 4 Dx Dy A_xy, so shifting A_xy(r) by 1/2 (pi/2) for all t is an
// exact symmetry. A slip of 1/2 between two time slices is the cheapest instanton (finite action,
// charge m = Dx Dy delta). spiderwebZ4Steps returns round(2 [e_xy]) per site, nonzero on such a slip.

#include <cmath>
#include <vector>

#include "MathUtil.h"
#include "Spiderweb.h"

namespace mcf {

struct SpiderwebFluxes {
	std::vector<double> e_xx, e_xy, b;   // wrapped onto (-1, 1], one value per site
};

inline SpiderwebFluxes spiderwebFluxes(const Spiderweb& s)
{
	const PeriodicLattice& lat = s.getLattice();
	const int n = lat.size();
	auto A = [&](int x, int y, int t, int type) {
		return s.variable(lat.index(lat.wrap(0, x), lat.wrap(1, y), lat.wrap(2, t)) * 3 + type);
	};
	SpiderwebFluxes f{ std::vector<double>(n), std::vector<double>(n), std::vector<double>(n) };
	for (int site = 0; site < n; site++)
	{
		const int x = lat.coord(site, 0), y = lat.coord(site, 1), t = lat.coord(site, 2);
		const double exx = -A(x, y, t + 1, 1) + A(x, y, t, 1)
			+ 4 * A(x + 1, y + 1, t, 0) - 4 * A(x + 1, y, t, 0) - 4 * A(x, y + 1, t, 0) + 4 * A(x, y, t, 0);
		const double exy = -A(x, y, t + 1, 2) + A(x, y, t, 2)
			+ A(x + 2, y, t, 0) - 2 * A(x + 1, y, t, 0) - A(x, y + 2, t, 0) + 2 * A(x, y + 1, t, 0);
		const double b = A(x + 2, y, t, 1) - 2 * A(x + 1, y, t, 1) - A(x, y + 2, t, 1) + 2 * A(x, y + 1, t, 1)
			- 4 * A(x + 1, y + 1, t, 2) + 4 * A(x + 1, y, t, 2) + 4 * A(x, y + 1, t, 2) - 4 * A(x, y, t, 2);
		f.e_xx[site] = mapToCircle(exx);
		f.e_xy[site] = mapToCircle(exy);
		f.b[site] = mapToCircle(b);
	}
	return f;
}

// Integer instanton density m(r), one value per site; the sum over the lattice is always 0.
inline std::vector<int> spiderwebInstantons(const Spiderweb& s)
{
	const PeriodicLattice& lat = s.getLattice();
	const SpiderwebFluxes f = spiderwebFluxes(s);
	auto at = [&](const std::vector<double>& v, int x, int y, int t) {
		return v[lat.index(lat.wrap(0, x), lat.wrap(1, y), lat.wrap(2, t))];
	};
	std::vector<int> m(lat.size());
	for (int site = 0; site < lat.size(); site++)
	{
		const int x = lat.coord(site, 0), y = lat.coord(site, 1), t = lat.coord(site, 2);
		const double twice_m = at(f.b, x, y, t + 1) - at(f.b, x, y, t)
			+ at(f.e_xx, x + 2, y, t) - 2 * at(f.e_xx, x + 1, y, t) - at(f.e_xx, x, y + 2, t) + 2 * at(f.e_xx, x, y + 1, t)
			- 4 * (at(f.e_xy, x + 1, y + 1, t) - at(f.e_xy, x + 1, y, t) - at(f.e_xy, x, y + 1, t) + at(f.e_xy, x, y, t));
		m[site] = (int)std::lround(twice_m / 2.0);
	}
	return m;
}

// Z4 steps: round(2 [e_xy]) per site, in {-2, ..., 2}; +-1 marks a pi/2 slip of A_xy at that site
// between slices t and t+1. An isolated relaxed step has |[e_xy]| ~ 0.32 at its core (KU = 0.5).
inline std::vector<int> spiderwebZ4Steps(const Spiderweb& s)
{
	const SpiderwebFluxes f = spiderwebFluxes(s);
	std::vector<int> steps(f.e_xy.size());
	for (size_t i = 0; i < steps.size(); i++)
		steps[i] = (int)std::lround(2.0 * f.e_xy[i]);
	return steps;
}

} // namespace mcf
