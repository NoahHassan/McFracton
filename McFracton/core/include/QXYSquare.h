#pragma once

#include <random>

#include "Lattice.h"
#include "System.h"

// Axes: x, y, tau. One angle per site.
//
// Out of scope for the architecture cleanup: it keeps compiling and is excluded from the
// regression tests. measure() reports the energy only.
class QXYSquare : public System {
public:
	QXYSquare(int size, int Ntau);
	QXYSquare(int size, int Ntau, float K_s, float K_t);
	~QXYSquare() override = default;
public:
	void randomize(std::mt19937& rng) override;
	double getEnergy() const override;
	double proposeUpdate(int index, double delta) const override;
	std::vector<std::pair<std::vector<int>, int>> getSpacialVortices() const;
	std::vector<std::pair<std::vector<int>, int>> getTemporalVortices() const;
	std::vector<std::string> observableNames() const override;
	std::vector<double> measure(double temperature) const override;
	std::vector<Channel> channels() const override;
	void fillChannel(int channel, std::vector<double>& out) const override;
public:
	const int size;
	const int Ntau;
	const int ss_size; // size of a space-space slice
	const int st_size; // size of a space-time slice
	const int nPlaqs;  // xy, x-tau and y-tau on every site
	float K_s;
	float K_t;
private:
	// first is real-space, second is tau-space
	std::pair<std::vector<int>, std::vector<int>> getSiteConnectedCluster(int siteIndex) const;
	std::pair<std::vector<int>, std::vector<int>> getPlaqConnectedCluster(int plaqIndex) const;
};
