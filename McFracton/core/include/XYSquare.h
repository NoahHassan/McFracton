#pragma once

#include "Lattice.h"
#include "System.h"

// Axes: x, y. One angle per site.
class XYSquare : public System {
public:
	XYSquare(int size);
	~XYSquare() override = default;
public:
	double getEnergy() const override;
	double getSinSqrX() const;
	double proposeUpdate(int index, double delta) const override;
	std::vector<std::pair<std::vector<int>, int>> getVortices() const;
	std::vector<std::string> observableNames() const override;
	std::vector<double> measure(double temperature) const override;
	std::vector<Channel> channels() const override;
	void fillChannel(int channel, std::vector<double>& out) const override;
public:
	const int size;
private:
	std::pair<std::vector<int>, std::vector<int>> getSiteConnectedCluster(int siteIndex) const;
	std::pair<std::vector<int>, std::vector<int>> getPlaqConnectedCluster(int plaqIndex) const;
};
