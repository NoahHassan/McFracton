#pragma once

#include <algorithm>

#include "Lattice.h"
#include "System.h"

class AbelianGaugeCube : public System {
public:
	AbelianGaugeCube(int linear_size, int temporal_size, unsigned int seed = 0);
	~AbelianGaugeCube() override = default;
public:
	double getEnergy() const;
	std::vector<double> getLocalEnergies() const { return {}; };
	double proposeSiteFlip(int index, double angle) const;
	double proposePlaqFlip(int index, double angle) const;
	void UpdateSite(int index, double angle);
	void UpdatePlaq(int index, double angle);
	virtual void OverrelaxSite(int index) override;
	std::vector<double> getFluxes_z() const;
	Observables Measure(double T) const;
public:
	const int linear_size;
	const int temporal_size;
	const int nSites;
	const int nPlaqs;
private:
	double getLocalEnergy_x(int nx, int ny, int nz, int nt, double angle) const;
	double getLocalEnergy_y(int nx, int ny, int nz, int nt, double angle) const;
	double getLocalEnergy_z(int nx, int ny, int nz, int nt, double angle) const;
	double getLocalEnergy_t(int nx, int ny, int nz, int nt, double angle) const;
	std::vector<std::pair<int, int>> getPlaqConnectedFields(int nx, int ny, int nz, int nt, int type) const;
	std::vector<std::pair<int, int>> getPlaqConnectedFields(int plaqIndex) const;
	double sum_plaquette(const std::vector<std::pair<int, int>>& plaquette) const;
	double sum_plaquette(const std::vector<std::pair<int, int>>& plaquette, int angle_index, double angle) const;
	int to_site_index(int nx, int ny, int nz, int nt) const;
	double get_field(int site_index, int direction) const;
	double get_field(int nx, int ny, int nz, int nt, int direction) const;
private:
	// Axes: x, y, z, t
	mcf::PeriodicLattice lattice;
	std::mt19937 rng;
	std::uniform_real_distribution<double> overrelax_dst;
};