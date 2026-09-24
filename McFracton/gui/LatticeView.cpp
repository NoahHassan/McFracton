#include "LatticeView.h"

#include <imgui.h>

#include <algorithm>
#include <cmath>

#include "ColorMaps.h"

void LatticeView::setSystem(const System& system)
{
	const mcf::PeriodicLattice& lattice = system.getLattice();

	axis_x = 0;
	axis_y = lattice.rank() > 1 ? 1 : 0;
	channel = 0;
	overlay = -1;
	slice.assign(lattice.rank(), 0);
	magnitude_scale = 1e-9;
	signed_scale = 1e-9;
	values.clear();
	overlay_values.clear();
}

void LatticeView::setArea(float x, float y, float width, float height)
{
	area_x = x;
	area_y = y;
	area_width = width;
	area_height = height;
}

std::vector<int> LatticeView::slicedAxes(const System& system) const
{
	std::vector<int> axes;
	for (int axis = 0; axis < system.getLattice().rank(); axis++)
		if (axis != axis_x && axis != axis_y)
			axes.push_back(axis);
	return axes;
}

void LatticeView::stepSlice(const System& system, int which, int offset)
{
	const std::vector<int> axes = slicedAxes(system);
	if (which < 0 || which >= (int)axes.size())
		return;

	const int axis = axes[which];
	const int extent = system.getLattice().extent(axis);
	slice[axis] = (slice[axis] + offset % extent + extent) % extent;
}

sf::Color LatticeView::colorOf(ChannelKind kind, double value) const
{
	switch (kind)
	{
	case ChannelKind::Angle:
		return mcf::NormalMapYellow(value);
	case ChannelKind::Integer:
		return mcf::RedWhiteBlue((int)std::lround(value));
	case ChannelKind::Signed:
		return mcf::RedWhiteBlue(value, signed_scale);
	case ChannelKind::Magnitude:
	default:
		return mcf::BlackWhite(value / magnitude_scale);
	}
}

void LatticeView::drawControls(const System& system)
{
	const mcf::PeriodicLattice& lattice = system.getLattice();
	const std::vector<Channel> channels = system.channels();
	if (channels.empty())
		return;

	// Channel, and the optional overlay drawn on top of it.
	std::vector<const char*> channel_names;
	for (const Channel& c : channels)
		channel_names.push_back(c.name.c_str());

	if (ImGui::Combo("Channel", &channel, channel_names.data(), (int)channel_names.size()))
	{
		magnitude_scale = 1e-9;
		signed_scale = 1e-9;
	}

	std::vector<const char*> overlay_names;
	overlay_names.push_back("none");
	for (const Channel& c : channels)
		overlay_names.push_back(c.name.c_str());
	int overlay_item = overlay + 1;
	if (ImGui::Combo("Overlay", &overlay_item, overlay_names.data(), (int)overlay_names.size()))
		overlay = overlay_item - 1;

	// Which plane of the lattice is on screen.
	if (lattice.rank() > 2)
	{
		std::vector<const char*> axis_names;
		for (const std::string& name : lattice.getAxisNames())
			axis_names.push_back(name.c_str());

		const int previous_x = axis_x;
		const int previous_y = axis_y;
		ImGui::Combo("Horizontal", &axis_x, axis_names.data(), (int)axis_names.size());
		ImGui::Combo("Vertical", &axis_y, axis_names.data(), (int)axis_names.size());
		// Two axes on the same screen direction would collapse the picture, so the one that did
		// not just change gives way.
		if (axis_x == axis_y)
		{
			if (axis_x != previous_x)
				axis_y = previous_x;
			else
				axis_x = previous_y;
		}
	}

	ImGui::SliderFloat("Zoom", &zoom, 0.3f, 1.25f, "%.2f");

	const std::vector<int> axes = slicedAxes(system);
	for (size_t i = 0; i < axes.size(); i++)
	{
		const int axis = axes[i];
		const std::string label = lattice.getAxisNames()[axis] + " slice";
		ImGui::SliderInt(label.c_str(), &slice[axis], 0, lattice.extent(axis) - 1);
	}
	if (!axes.empty())
		ImGui::TextDisabled("left/right step %s, enter steps %s",
			lattice.getAxisNames()[axes[0]].c_str(),
			axes.size() > 1 ? lattice.getAxisNames()[axes[1]].c_str() : "nothing");
}

void LatticeView::draw(const System& system)
{
	const mcf::PeriodicLattice& lattice = system.getLattice();
	const std::vector<Channel> channels = system.channels();
	if (channels.empty() || channel < 0 || channel >= (int)channels.size())
		return;

	system.fillChannel(channel, values);
	if ((int)values.size() != lattice.size())
		return;

	const bool has_overlay = overlay >= 0 && overlay < (int)channels.size();
	if (has_overlay)
	{
		system.fillChannel(overlay, overlay_values);
		if ((int)overlay_values.size() != lattice.size())
			return;
	}

	// The scales for the two kinds that have no fixed range grow to fit what has been seen.
	const ChannelKind kind = channels[channel].kind;
	if (kind == ChannelKind::Magnitude || kind == ChannelKind::Signed)
	{
		double largest = 0.0;
		for (double v : values)
			largest = std::max(largest, std::abs(v));
		if (kind == ChannelKind::Magnitude)
			magnitude_scale = std::max(magnitude_scale, largest);
		else
			signed_scale = std::max(signed_scale, largest);
	}

	const int nx = lattice.extent(axis_x);
	const int ny = lattice.extent(axis_y);
	if (nx <= 0 || ny <= 0)
		return;

	// Square cells, centred in the area. The lattice takes a fraction of the space rather than
	// all of it, so it sits in the window instead of filling it edge to edge; zoom scales that.
	const float fill = 0.78f * zoom;
	const float cell = std::max(2.0f, std::min(area_width * fill / nx, area_height * fill / ny));
	const float origin_x = area_x + (area_width - cell * nx) / 2.0f;
	const float origin_y = area_y + (area_height - cell * ny) / 2.0f;
	// The grid line between cells, kept proportional so it reads the same at any lattice size,
	// and dropped once the cells are too small to spare the pixels.
	const float gap = cell > 6.0f ? std::clamp(cell * 0.10f, 2.0f, 8.0f) : 0.0f;

	cells.setPrimitiveType(sf::Quads);
	cells.resize((size_t)nx * ny * 4);

	coords = slice;
	size_t vertex = 0;
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			coords[axis_x] = i;
			coords[axis_y] = j;
			int site = 0;
			for (int axis = 0; axis < lattice.rank(); axis++)
				site += coords[axis] * lattice.stride(axis);

			sf::Color color = colorOf(kind, values[site]);
			if (has_overlay)
			{
				// An overlay only marks the cells where it is nonzero, leaving the rest of the
				// picture alone. This is what "draw the monopoles on top" used to do by hand.
				const double marker = overlay_values[site];
				if (marker != 0.0)
					color = mcf::RedWhiteBlue(marker > 0.0 ? 1 : -1);
			}

			const float x0 = origin_x + i * cell;
			const float y0 = origin_y + j * cell;
			const float x1 = x0 + cell - gap;
			const float y1 = y0 + cell - gap;

			cells[vertex + 0] = sf::Vertex({ x0, y0 }, color);
			cells[vertex + 1] = sf::Vertex({ x1, y0 }, color);
			cells[vertex + 2] = sf::Vertex({ x1, y1 }, color);
			cells[vertex + 3] = sf::Vertex({ x0, y1 }, color);
			vertex += 4;
		}
	}

	window.draw(cells);
}
