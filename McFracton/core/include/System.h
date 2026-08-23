#pragma once

#include <random>

class System {
public:
	struct Observables {
		double energy = 0.0;
		double helicity_modulus = 0.0;
		double polyakov_loop = 0.0;
		double flux_cos = 0.0;
		int n_defects_a = 0;
		int n_defects_b = 0;
	};
public:
	System(int n_site_variables, int n_plaq_variables)
		: 
		n_site_variables(n_site_variables), n_plaq_variables(n_plaq_variables)
	{};
	virtual ~System() = 0;
public:
	virtual double getEnergy() const = 0;
	virtual double proposeSiteFlip(int index, double angle) const = 0;
	virtual double proposePlaqFlip(int index, double angle) const = 0;
	virtual void UpdateSite(int index, double angle) = 0;
	virtual void UpdatePlaq(int index, double angle) = 0;
	virtual void OverrelaxSite(int index) { throw("OverrelaxSite not implemented"); };
	virtual void OverrelaxPlaq(int index) { throw("OverrelaxPlaq not implemented"); };
	double getSite(int index) const;
	double getPlaq(int index) const;
	virtual Observables Measure(double T) const { throw("Measure not implemented"); return {}; };
public:
	const int n_site_variables;
	const int n_plaq_variables;
protected:
	std::vector<double> site_fields;
	std::vector<double> plaq_fields;
};