#pragma once

#include <SFML/Graphics.hpp>
#include <assert.h>

#include "XYSquare.h"
#include "QXYSquare.h"
#include "AbelianGaugeSquare.h"
#include "Square.h"

#define PI 3.1415926535897932384

class Canvas {
public:
	Canvas(sf::RenderWindow& window)
		:
		window(window)
	{}
	Canvas(sf::RenderWindow& window, Vec2D offset)
		:
		window(window), offset(offset)
	{}
	void Initialize(const AbelianGaugeSquare& field, float square_size)
	{
		site_pixels = std::vector<Square>(field.linear_size * field.linear_size);
		for (int n = 0; n < field.linear_size * field.linear_size; n++)
		{
			int ny = n / field.linear_size;
			int nx = n - ny * field.linear_size;

			Square sq(
				Vec2D(
					square_size * (nx - field.linear_size / 2) + window.getSize().x / 2,
					square_size * (ny - field.linear_size / 2) + window.getSize().y / 2
				) + offset,
				square_size
			);

			sq.SetOutlineThickness(2.0);
			sq.SetOutlineColor(sf::Color::Black);

			site_pixels[n] = std::move(sq);
		}
	}
	void Draw(const AbelianGaugeSquare& field, int direction)
	{
		UpdateFieldColors(field, direction);
		for (const auto& sq : site_pixels)
		{
			window.draw(sq);
		}
	}
	void Draw(const AbelianGaugeSquare& field, int direction, int layer)
	{
		UpdateFieldColors(field, direction, layer);
		for (const auto& sq : site_pixels)
		{
			window.draw(sq);
		}
	}
	void DrawMonopoles(const AbelianGaugeSquare& field, int layer)
	{
		//Interpret pixels as plaquettes now (nSites = nPlaquettes)
		UpdateMonopoleColors(field, layer);
		for (const auto& sq : site_pixels)
		{
			window.draw(sq);
		}
	}
	void DrawFluxes(const AbelianGaugeSquare& field, int layer)
	{
		UpdateFluxColors(field, layer);
		for (const auto& sq : site_pixels)
		{
			window.draw(sq);
		}
	}
private:
	void UpdateFieldColors(const AbelianGaugeSquare& field, int direction)
	{
		for (int n = 0; n < field.linear_size * field.linear_size; n++)
		{
			int field_index = n * 3 + direction;
			assert(field_index < field.numVariables());
			const double theta = field.variable(field_index);

			//site_pixels[n].SetFillColor(NormalMapYellow(theta, 0.3, 0.6, 0.8));
			site_pixels[n].SetFillColor(GreenRedUniform(theta));
		}
	}
	void UpdateFieldColors(const AbelianGaugeSquare& field, int direction, int layer)
	{
		for (int n = 0; n < field.linear_size * field.linear_size; n++)
		{
			int n_shifted = n + (int)site_pixels.size() * layer;
			int field_index = n_shifted * 3 + direction;
			assert(field_index < field.numVariables());
			const double theta = field.variable(field_index);

			site_pixels[n].SetFillColor(NormalMapYellow(theta, 0.3, 0.6, 0.8));
			//site_pixels[n].SetFillColor(GreenRedUniform(theta));
		}
	}
	void UpdateMonopoleColors(const AbelianGaugeSquare& field, int layer)
	{
		const std::vector<int> monopoles = field.getMonopoles();
		for (int n = 0; n < field.linear_size * field.linear_size; n++)
		{
			int n_shifted = n + (int)site_pixels.size() * layer;
			if(monopoles[n_shifted] != 0.0)
				site_pixels[n].SetFillColor(RedWhiteBlue(monopoles[n_shifted]));
		}
	}
	void UpdateFluxColors(const AbelianGaugeSquare& field, int layer)
	{
		const std::vector<double> fluxes = field.getFluxes_z();
		for (int n = 0; n < field.linear_size * field.linear_size; n++)
		{
			int n_shifted = n + (int)site_pixels.size() * layer;
			//site_pixels[n].SetFillColor(NormalMapYellow(fluxes[n_shifted], 0.3, 0.6, 0.8));
			site_pixels[n].SetFillColor(BlackWhite(fluxes[n_shifted]));
		}
	}
	sf::Color GreenRedUniform(const double& theta)
	{
		const int theta_c = int(std::abs(theta) * 510.0) % 510;
		sf::Uint8 r = sf::Uint8(255 - std::abs(theta_c - 255));
		sf::Uint8 g = sf::Uint8(std::abs(theta_c - 255));
		return sf::Color(r, g, 122u);
	}
	sf::Color BlackWhite(const double& cosine)
	{
		max_flux = std::max(abs(cosine), max_flux);

		const double mapped_cos = abs(cosine / max_flux);
		const int color_val = int(mapped_cos * 255) % 255;
		return sf::Color(color_val, color_val, color_val);

	}
	sf::Color NormalMapYellow(const double& theta, double tilt, double intensity, double ambient)
	{
		double lightAngle = intensity * (1.0 + cos(2.0 * PI * (theta - tilt))) / 0.5;

		sf::Uint8 r = sf::Uint8(std::min(255.0, (lightAngle + ambient) * 90));
		sf::Uint8 g = sf::Uint8(std::min(255.0, (lightAngle + ambient) * 60));
		sf::Uint8 b = sf::Uint8(std::min(255.0, (lightAngle + ambient) * 20));

		return sf::Color(r, g, b);
	}
	sf::Color RedWhiteBlue(int n)
	{
		if (n > 0)
			return sf::Color::Blue;
		if (n < 0)
			return sf::Color::Red;
		return sf::Color::White;
	}
private:
	sf::RenderWindow& window;
	std::vector<Square> site_pixels;
	double max_flux = 0.0;
	Vec2D offset{ 0.0f, 0.0f };
};