#pragma once

#include <SFML/Graphics.hpp>

#include <string>
#include <vector>

#include "System.h"

/// <summary>
/// The one canvas, replacing Canvas, HyperCanvas and SpiderCanvas.
///
/// It knows nothing about any particular system. It asks the system for its lattice, its channels
/// and the values of one channel, then draws a plane through that lattice: two axes are shown, and
/// every remaining axis is held at a slice the user can move. That covers 2D, 2+1D and 3+1D alike,
/// so a new system needs no GUI work at all beyond registering itself.
///
/// Cells are one sf::VertexArray of quads with a one pixel gap, which is what gives the grid look
/// without a shape object per site.
/// </summary>
class LatticeView {
public:
	explicit LatticeView(sf::RenderWindow& window) : window(window) {}

	/// <summary>
	/// Point the view at a system. Resets the axes, the channel and every slice, so it is also what
	/// is called after the user rebuilds with different parameters.
	/// </summary>
	void setSystem(const System& system);

	/// <summary>
	/// The ImGui block for what is on screen: the channel, the overlay, which two axes are drawn
	/// and where the remaining ones are sliced.
	/// </summary>
	void drawControls(const System& system);

	void draw(const System& system);

	/// <summary>
	/// Steps one of the sliced axes, for the arrow keys and Enter. `which` is 0 for the first
	/// sliced axis and 1 for the second; out of range or absent axes are ignored.
	/// </summary>
	void stepSlice(const System& system, int which, int offset);

	/// <summary>
	/// The screen rectangle the lattice is drawn into. Main sets this to whatever the control
	/// panel leaves over.
	/// </summary>
	void setArea(float x, float y, float width, float height);

private:
	// The axes that are not displayed, in ascending order. These are the ones with slice sliders.
	std::vector<int> slicedAxes(const System& system) const;
	// Colour for one value of the current channel.
	sf::Color colorOf(ChannelKind kind, double value) const;

private:
	sf::RenderWindow& window;
	sf::VertexArray cells{ sf::Quads };

	std::vector<double> values;
	std::vector<double> overlay_values;
	std::vector<int> coords;
	// One entry per axis; only the entries for non-displayed axes are used.
	std::vector<int> slice;

	int axis_x = 0;
	int axis_y = 1;
	int channel = 0;
	int overlay = -1;   // -1 means no overlay.
	// Scales how much of the available area the lattice takes, so the default size is a starting
	// point rather than something to recompile.
	float zoom = 1.0f;

	// Magnitude channels have no natural range, so the scale grows to the largest value seen and
	// is reset whenever the channel or the system changes - the running maximum the old canvas kept.
	double magnitude_scale = 1e-9;
	double signed_scale = 1e-9;

	float area_x = 0.0f;
	float area_y = 0.0f;
	float area_width = 0.0f;
	float area_height = 0.0f;
};
