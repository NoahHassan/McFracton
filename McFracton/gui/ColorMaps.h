#pragma once

#include <SFML/Graphics.hpp>

#include <algorithm>
#include <cmath>

#include "MathUtil.h"

/// <summary>
/// The colour maps the three old canvases each carried their own identical copy of, moved here
/// unchanged. A channel's ChannelKind picks which one the view uses, so a system never says
/// anything about colour.
/// </summary>
namespace mcf {

	/// <summary>
	/// Green to red and back, repeating. Period 2 in theta, matching the field convention.
	/// </summary>
	inline sf::Color GreenRedUniform(double theta)
	{
		const int theta_c = int(std::abs(theta) * 510.0) % 510;
		sf::Uint8 r = sf::Uint8(255 - std::abs(theta_c - 255));
		sf::Uint8 g = sf::Uint8(std::abs(theta_c - 255));
		return sf::Color(r, g, 122u);
	}

	/// <summary>
	/// Black for the lowest magnitudes, white for the highest. The caller divides by its own scale.
	/// </summary>
	inline sf::Color BlackWhite(double value)
	{
		// map cos in range (-3,3) to range (0,2) and apply triangle function
		const double mapped_cos = std::abs(value / 2.0 + 0.5); // lowest energies colored black, highest white
		const int color_val = int(mapped_cos * 510) % 510;
		return sf::Color(sf::Uint8(color_val), sf::Uint8(color_val), sf::Uint8(color_val));
	}

	/// <summary>
	/// An angle lit as if it were a surface normal, which is what makes vortices visible as
	/// pinwheels rather than as colour noise.
	/// </summary>
	inline sf::Color NormalMapYellow(double theta, double tilt = 0.3, double intensity = 0.6, double ambient = 0.8)
	{
		double lightAngle = intensity * (1.0 + std::cos(2.0 * kPi * (theta - tilt))) / 0.5;

		sf::Uint8 r = sf::Uint8(std::min(255.0, (lightAngle + ambient) * 90));
		sf::Uint8 g = sf::Uint8(std::min(255.0, (lightAngle + ambient) * 60));
		sf::Uint8 b = sf::Uint8(std::min(255.0, (lightAngle + ambient) * 20));

		return sf::Color(r, g, b);
	}

	/// <summary>
	/// Sign only: blue positive, red negative, white zero. For genuinely integer channels such as
	/// vortices and monopoles, where the count is small and the sign is the whole story.
	/// </summary>
	inline sf::Color RedWhiteBlue(int n)
	{
		if (n > 0)
			return sf::Color::Blue;
		if (n < 0)
			return sf::Color::Red;
		return sf::Color::White;
	}

	/// <summary>
	/// The continuous version of the above, for signed channels like the flux, where the magnitude
	/// carries as much as the sign and three flat colours would throw it away. White at zero,
	/// saturating to the same blue and red at +-scale.
	/// </summary>
	inline sf::Color RedWhiteBlue(double value, double scale)
	{
		if (scale <= 0.0)
			return sf::Color::White;

		const double t = std::clamp(value / scale, -1.0, 1.0);
		const double magnitude = std::abs(t);
		// Interpolate white -> blue for positive, white -> red for negative.
		const sf::Uint8 faded = sf::Uint8(255.0 * (1.0 - magnitude));
		if (t >= 0.0)
			return sf::Color(faded, faded, 255u);
		return sf::Color(255u, faded, faded);
	}

} // namespace mcf
