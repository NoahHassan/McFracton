#pragma once

#include <algorithm>
#include <array>
#include <random>

#include "Lattice.h"
#include "System.h"

// Axes: x, y, t. Rank-2 fields on every site, so the variable index is site * 3 + type:
// 0: A_0, 1: A_xx, 2: A_xy.
//
// Overrelaxation is deliberately not implemented, so it falls back to the base class throw.
class Spiderweb : public System {
public:
	Spiderweb(int linear_size, int temporal_size, double KU, unsigned int seed = 0);
	~Spiderweb() override = default;
public:
	double getEnergy() const override;
	std::vector<double> getLocalEnergies() const;
	double getEnergy(std::vector<double>& localFluxes) const;
	double proposeUpdate(int index, double delta) const override;
	std::vector<std::string> observableNames() const override;
	std::vector<double> measure(double temperature) const override;
	std::vector<Channel> channels() const override;
	void fillChannel(int channel, std::vector<double>& out) const override;
public:
	const int linear_size;
	const int spatial_size;
	const int temporal_size;
	const int nSites;
	const double KU;
private:
	// Fills one local energy per site and returns their sum; the single implementation behind
	// getEnergy(), getLocalEnergies() and getEnergy(std::vector<double>&).
	double accumulateLocalEnergies(std::vector<double>& localEnergies) const;
	std::vector<std::pair<int, double>> getElectricTerms_xx(int site_index) const;
	std::vector<std::pair<int, double>> getElectricTerms_xy(int site_index) const;
	std::vector<std::pair<int, double>> getMagneticTerms(int site_index) const;
	std::array<int, 3> index_from_site(int site_index) const;
	int field_index_from_site(int nx, int ny, int nt, int type) const;
	int field_index_from_site(int site_index, int type) const;
	double get_field(int site_index, int type) const;
	double get_field(int nx, int ny, int nt, int type) const;
private:
	std::mt19937 rng;
	std::uniform_real_distribution<double> overrelax_dst;
};
