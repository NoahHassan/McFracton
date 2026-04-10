#pragma once

#include <fstream>

#include "System.h"

class QXYSquare : public System {
public:
	QXYSquare(int size, int Ntau);
	QXYSquare(int size, int Ntau, float K_s, float K_t);
	~QXYSquare() override = default;
public:
	double getEnergy() const;
	double proposeSiteFlip(int index, double angle) const;
	double proposePlaqFlip(int index, double angle) const;
	void UpdateSite(int index, double angle);
	void UpdatePlaq(int index, double angle);
	std::vector<std::pair<std::vector<int>, int>> getVortices() const;
	void LogToFile(std::ofstream& outfile) const;
public:
	const int size;
	const int Ntau;
	float K_s;
	float K_t;
private:
	// first is real-space, second is tau-space
	const std::pair<std::vector<int>, std::vector<int>> getSiteConnectedCluster(int siteIndex) const;
	double mapToCircle(const double& d) const;

	int sgn(double val) const {
		return (0.0 < val) - (val < 0.0);
	}
};