#pragma once

#include <array>
#include <cassert>
#include <string>
#include <vector>

namespace mcf {

	/// <summary>
	/// A rectangular lattice with periodic boundaries in every direction.
	///
	/// Sites are numbered row-major with axis 0 running fastest, which is the layout every system
	/// already used by hand: index = nx + L*(ny + L*(nz + ...)). The class only replaces the index
	/// arithmetic and the "(n - 1 + L) % L" wraps that were copied into each system; it holds no
	/// field data and computes exactly the same integers as the code it replaces.
	/// </summary>
	class PeriodicLattice {
	public:
		PeriodicLattice() = default;

		PeriodicLattice(std::vector<int> extents, std::vector<std::string> axis_names)
			:
			extents(std::move(extents)), axis_names(std::move(axis_names))
		{
			assert(this->extents.size() == this->axis_names.size());
			strides.resize(this->extents.size());
			int stride = 1;
			for (size_t axis = 0; axis < this->extents.size(); axis++)
			{
				assert(this->extents[axis] > 0);
				strides[axis] = stride;
				stride *= this->extents[axis];
			}
			total_sites = stride;
		}

	public:
		int rank() const { return (int)extents.size(); }
		int extent(int axis) const { return extents[axis]; }
		int stride(int axis) const { return strides[axis]; }
		int size() const { return total_sites; }
		const std::vector<int>& getExtents() const { return extents; }
		const std::vector<std::string>& getAxisNames() const { return axis_names; }

		/// <summary>
		/// Site index from coordinates, axis 0 first: index(nx, ny, nt). Coordinates must already
		/// be in range; use wrap() for neighbours that step off the edge.
		/// </summary>
		template<typename... Coords>
		int index(Coords... coords) const
		{
			const std::array<int, sizeof...(Coords)> c{ coords... };
			assert((int)c.size() == rank());
			int result = 0;
			for (size_t axis = 0; axis < c.size(); axis++)
			{
				assert(0 <= c[axis] && c[axis] < extents[axis]);
				result += c[axis] * strides[axis];
			}
			return result;
		}

		/// <summary>
		/// Folds a coordinate back into [0, extent) for a step of at most one lattice length,
		/// which is what "(n + 1) % L" and "(n - 1 + L) % L" did.
		/// </summary>
		int wrap(int axis, int coord) const
		{
			const int L = extents[axis];
			assert(-L <= coord && coord < 2 * L);
			return (coord + L) % L;
		}

		/// <summary>
		/// The neighbour of a site along one axis, offset steps away (offset may be negative).
		/// </summary>
		int neighbor(int site_index, int axis, int offset) const
		{
			const int L = extents[axis];
			const int old_coord = coord(site_index, axis);
			const int new_coord = wrap(axis, old_coord + offset);
			return site_index + (new_coord - old_coord) * strides[axis];
		}

		/// <summary>
		/// One coordinate of a site index.
		/// </summary>
		int coord(int site_index, int axis) const
		{
			return (site_index / strides[axis]) % extents[axis];
		}

		void coords(int site_index, std::vector<int>& out) const
		{
			out.resize(extents.size());
			for (size_t axis = 0; axis < extents.size(); axis++)
				out[axis] = coord(site_index, (int)axis);
		}

	private:
		std::vector<int> extents;
		std::vector<std::string> axis_names;
		std::vector<int> strides;
		int total_sites = 0;
	};

} // namespace mcf
