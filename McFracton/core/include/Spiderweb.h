#pragma once

#include <fstream>
#include <algorithm>
#include <array>

#include "System.h"

// On each site are attached:
// 1 link in x-direction
// 1 link in y-direction
// 1 link in t-direction
// 1 xy-plaquette

class Spiderweb : public System {
public:
	Spiderweb(int linear_size, int temporal_size, double KU);
	~Spiderweb() override = default;
public:
	double getEnergy() const;
	double proposeSiteFlip(int index, double angle) const;
	double proposePlaqFlip(int index, double angle) const;
	void UpdateSite(int index, double angle);
	void UpdatePlaq(int index, double angle);
	virtual void OverrelaxSite(int index) override;
	Observables Measure(double T) const;
public:
	const int linear_size;
	const int spatial_size;
	const int temporal_size;
	const int nSites;
	const int nPlaqs;
	const double KU;
private:
	std::vector<std::pair<int, double>> getElectricTerms_xx(int site_index) const;
	std::vector<std::pair<int, double>> getElectricTerms_xy(int site_index) const;
	std::vector<std::pair<int, double>> getMagneticTerms(int site_index) const;
	const int to_site_index(int nx, int ny, int nt) const;
	const std::array<int, 3> index_from_site(int site_index) const;
	const int field_index_from_site(int nx, int ny, int nt, int type) const;
	const int field_index_from_site(int site_index, int type) const;
	double get_field(int site_index, int type) const;
	double get_field(int nx, int ny, int nt, int type) const;
	double mapToCircle(const double& d) const;

	int sgn(double val) const {
		return (0.0 < val) - (val < 0.0);
	}
private:
	std::mt19937 rng;
	std::uniform_real_distribution<double> overrelax_dst;
};