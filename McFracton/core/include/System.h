#pragma once

#include <random>

class System {
public:
	System(int nSites, int nPlaqs) : nSites(nSites), nPlaqs(nPlaqs) {};
	virtual ~System() = 0;
public:
	virtual double getEnergy() const = 0;
	virtual double proposeSiteFlip(int index, double angle) const = 0;
	virtual double proposePlaqFlip(int index, double angle) const = 0;
	virtual void UpdateSite(int index, double angle) = 0;
	virtual void UpdatePlaq(int index, double angle) = 0;
	double getSite(int index) const;
	double getPlaq(int index) const;
	virtual void LogToFile(std::ofstream& outfile) const = 0;
protected:
	virtual const std::pair<std::vector<int>, std::vector<int>> getSiteConnectedCluster(int siteIndex) const = 0;
	virtual const std::pair<std::vector<int>, std::vector<int>> getPlaqConnectedCluster(int plaqIndex) const = 0;
public:
	const int nSites;
	const int nPlaqs;
protected:
	std::vector<double> site_fields;
	std::vector<double> plaq_fields;
};