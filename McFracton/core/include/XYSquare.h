#pragma once

#include "Lattice.h"
#include "System.h"

class XYSquare : public System {
public:
	XYSquare(int size);
	~XYSquare() override = default;
public:
	double getEnergy() const;
	std::vector<double> getLocalEnergies() const { return {}; };
	double getSinSqrX() const;
	double proposeSiteFlip(int index, double angle) const;
	double proposePlaqFlip(int index, double angle) const;
	void UpdateSite(int index, double angle);
	void UpdatePlaq(int index, double angle);
	std::vector<std::pair<std::vector<int>, int>> getVortices() const;
	Observables Measure(double T) const;
public:
	const int size;
private:
	std::pair<std::vector<int>, std::vector<int>> getSiteConnectedCluster(int siteIndex) const;
	std::pair<std::vector<int>, std::vector<int>> getPlaqConnectedCluster(int plaqIndex) const;
private:
	// Axes: x, y
	mcf::PeriodicLattice lattice;
};