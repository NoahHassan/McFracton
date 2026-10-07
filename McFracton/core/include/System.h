#pragma once

#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "Lattice.h"

/// <summary>
/// What a visualisation channel means, which is all the GUI needs in order to pick a colormap.
/// </summary>
enum class ChannelKind { Angle, Signed, Integer, Magnitude };

/// <summary>
/// One scalar field over the lattice that a system can hand out for drawing.
///
/// low and high are the colour range of a Signed or Magnitude channel, for a system that knows it
/// in advance: a Signed channel is white at 0 and saturates at +-high, a Magnitude channel runs
/// from black at low to white at high. Left equal, the system declares no range, and the view
/// measures one from the state it is first shown and then holds it.
/// </summary>
struct Channel {
	std::string name;
	ChannelKind kind;
	double low = 0.0;
	double high = 0.0;
};

/// <summary>
/// A system is a set of real variables on a periodic lattice, plus the energy differences that
/// Metropolis needs, the observables it wants logged, and the fields it wants drawn.
///
/// Observables are declared by name once (observableNames) and returned as plain doubles in the
/// same order (measure), so McMachine can build its log header and its statistics without knowing
/// which system it is running. Visualisation is pure data (channels / fillChannel), so the core
/// stays free of SFML.
///
/// A system owns no random number generator. Every random draw it needs comes from the generator
/// McMachine passes in, so one seed fixes the whole run: the initial configuration and the Markov
/// chain alike.
/// </summary>
class System {
public:
	System(mcf::PeriodicLattice site_lattice, int n_variables)
		:
		lattice(std::move(site_lattice)), fields(n_variables), n_variables(n_variables)
	{};
	virtual ~System() = 0;
public:
	int numVariables() const { return n_variables; }
	double variable(int index) const;
	const mcf::PeriodicLattice& getLattice() const { return lattice; }

	// Initial configuration. McMachine calls this once, before the first sweep. The default is a
	// cold start (every variable 0), which is what XYSquare has always started from.
	virtual void randomize(std::mt19937& rng) {};
	// Puts the system into a ground state, for runs that heat up from it. A system that has not
	// said what its ground state is has to fail loudly rather than start from something else.
	virtual void setGroundState() { throw std::logic_error("setGroundState not implemented"); };

	// Metropolis
	virtual double getEnergy() const = 0;
	virtual double proposeUpdate(int index, double delta) const = 0;
	void applyUpdate(int index, double delta);
	// Unimplemented overrelaxation has to fail loudly rather than silently do nothing.
	virtual void overrelax(int index, std::mt19937& rng) { throw std::logic_error("overrelax not implemented"); };

	// Measurement. The names are the log-file column names; measure() returns one value per name,
	// in the same order.
	virtual std::vector<std::string> observableNames() const = 0;
	virtual std::vector<double> measure(double temperature) const = 0;

	// Visualisation. fillChannel writes one value per lattice site, so out.size() == lattice.size().
	// A site the channel's field does not live on gets NaN, which the view leaves undrawn.
	virtual std::vector<Channel> channels() const = 0;
	virtual void fillChannel(int channel, std::vector<double>& out) const = 0;
protected:
	mcf::PeriodicLattice lattice;
	std::vector<double> fields;
private:
	const int n_variables;
};
