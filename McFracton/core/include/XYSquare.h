#pragma once

#include "System.h"

class XYSquare : public System {
public:
	XYSquare(int size);
	~XYSquare() override = default;
public:
	double getEnergy() const;
	double proposeSiteFlip(int index, double angle) const;
	double proposePlaqFlip(int index, double angle) const;
	void UpdateSite(int index, double angle);
	void UpdatePlaq(int index, double angle);
	std::vector<std::pair<std::vector<int>, int>> getVortices() const;
public:
	const int length;
private:
	const std::pair<std::vector<int>, std::vector<int>> getSiteConnectedCluster(int siteIndex) const;
	const std::pair<std::vector<int>, std::vector<int>> getPlaqConnectedCluster(int plaqIndex) const;
	double mapToCircle(const double& d) const;

	int sgn(double val) const {
		return (0.0 < val) - (val < 0.0);
	}
};