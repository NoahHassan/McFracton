#pragma once

#include <algorithm>
#include <random>

#include "Lattice.h"
#include "System.h"

// Axes: x, y, t. A_x, A_y, A_t on every site, so the variable index is site * 3 + direction.
class AbelianGaugeSquare : public System {
public:
	AbelianGaugeSquare(int linear_size, int temporal_size);
	~AbelianGaugeSquare() override = default;
public:
	void randomize(std::mt19937& rng) override;
	double getEnergy() const override;
	double proposeUpdate(int index, double delta) const override;
	void overrelax(int index, std::mt19937& rng) override;
	std::vector<int> getMonopoles() const;
	std::vector<double> getFluxes_z() const;
	std::vector<std::string> observableNames() const override;
	std::vector<double> measure(double temperature) const override;
	std::vector<Channel> channels() const override;
	void fillChannel(int channel, std::vector<double>& out) const override;
public:
	const int linear_size;
	const int temporal_size;
	const int nSites;
	const int nPlaqs;
private:
	double getLocalEnergy_x(int nx, int ny, int nt, double angle) const;
	double getLocalEnergy_y(int nx, int ny, int nt, double angle) const;
	double getLocalEnergy_t(int nx, int ny, int nt, double angle) const;
	std::vector<std::pair<int, int>> getPlaqConnectedFields(int nx, int ny, int nt, int type) const;
	std::vector<std::pair<int, int>> getPlaqConnectedFields(int plaqIndex) const;
	double sum_plaquette(const std::vector<std::pair<int, int>>& plaquette) const;
	double sum_plaquette(const std::vector<std::pair<int, int>>& plaquette, int angle_index, double angle) const;
	int to_site_index(int nx, int ny, int nt) const;
	double get_field(int site_index, int direction) const;
	double get_field(int nx, int ny, int nt, int direction) const;
};
