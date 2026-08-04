#pragma once

#include <fstream>
#include <algorithm>

#include "System.h"

class AbelianGaugeSquare : public System {
public:
	AbelianGaugeSquare(int linear_size, int temporal_size);
	~AbelianGaugeSquare() override = default;
public:
	double getEnergy() const;
	double proposeSiteFlip(int index, double angle) const;
	double proposePlaqFlip(int index, double angle) const;
	void UpdateSite(int index, double angle);
	void UpdatePlaq(int index, double angle);
	virtual void OverrelaxSite(int index) override;
	std::vector<int> getMonopoles() const;
	std::vector<double> getFluxes_z() const;
	Observables Measure(double T) const;
	//void LogToFile(std::ofstream& outfile) const;
public:
	const int linear_size;
	const int temporal_size;
	const int nSites;
	const int nPlaqs;
private:
	double getLocalEnergy_x(int nx, int ny, int nt, double angle) const;
	double getLocalEnergy_y(int nx, int ny, int nt, double angle) const;
	double getLocalEnergy_t(int nx, int ny, int nt, double angle) const;
	const std::vector<std::pair<int, int>> getPlaqConnectedFields(int nx, int ny, int nt, int type) const;
	const std::vector<std::pair<int, int>> getPlaqConnectedFields(int plaqIndex) const;
	const double sum_plaquette(const std::vector<std::pair<int, int>>& plaquette) const;
	const double sum_plaquette(const std::vector<std::pair<int, int>>& plaquette, int angle_index, double angle) const;
	int to_site_index(int nx, int ny, int nt) const;
	double get_field(int site_index, int direction) const;
	double get_field(int nx, int ny, int nt, int direction) const;
	double mapToCircle(const double& d) const;

	int sgn(double val) const {
		return (0.0 < val) - (val < 0.0);
	}
private:
	std::mt19937 rng;
	std::uniform_real_distribution<double> overrelax_dst;
};