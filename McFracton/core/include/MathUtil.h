#pragma once

#include <cmath>
#include <numbers>

namespace mcf {

	inline constexpr double kPi = std::numbers::pi;

	// Sign of a value: -1, 0 or +1.
	inline int sgn(double value)
	{
		return (0.0 < value) - (value < 0.0);
	}

	/// <summary>
	/// Wraps a value onto the interval (-1, 1], the period of the fields, which are measured in
	/// units where the action is cos(pi * x). Used to fold plaquette sums back before counting
	/// defects, so that a 2 pi jump does not register as one.
	/// </summary>
	inline double mapToCircle(double d)
	{
		const double half = d / 2.0;
		const double wrapped_half = (half >= 0.0) ? half - int(half + 0.5) : half - int(half - 0.5);
		return 2.0 * wrapped_half;
	}

}