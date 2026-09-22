#pragma once

#include <SFML/Graphics.hpp>
#include <assert.h>

#include "XYSquare.h"
#include "QXYSquare.h"
#include "AbelianGaugeCube.h"
#include "Square.h"

#define PI 3.1415926535897932384

class HyperCanvas {
public:
	HyperCanvas(sf::RenderWindow& window)
		:
		window(window)
	{}
	HyperCanvas(sf::RenderWindow& window, Vec2D offset)
		:
		window(window), offset(offset)
	{}
	void Initialize(const AbelianGaugeCube& field, float square_size)
	{
		planar_size = field.linear_size * field.linear_size;
		cubic_size = planar_size * field.linear_size;

		site_pixels = std::vector<Square>(planar_size);
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
	void Draw(const AbelianGaugeCube& field, int direction)
	{
		UpdateFieldColors(field, direction);
		for (const auto& sq : site_pixels)
		{
			window.draw(sq);
		}
	}
	void Draw(const AbelianGaugeCube& field, int direction, int layer, int time)
	{
		UpdateFieldColors(field, direction, layer, time);
		for (const auto& sq : site_pixels)
		{
			window.draw(sq);
		}
	}
	void DrawFluxes(const AbelianGaugeCube& field, int layer, int time)
	{
		UpdateFluxColors(field, layer, time);
		for (const auto& sq : site_pixels)
		{
			window.draw(sq);
		}
	}
private:
	void UpdateFieldColors(const AbelianGaugeCube& field, int direction)
	{
		for (int n = 0; n < field.linear_size * field.linear_size; n++)
		{
			int field_index = n * 4 + direction;
			assert(field_index < field.numVariables());
			const double theta = field.variable(field_index);

			//site_pixels[n].SetFillColor(NormalMapYellow(theta, 0.3, 0.6, 0.8));
			site_pixels[n].SetFillColor(GreenRedUniform(theta));
		}
	}
	void UpdateFieldColors(const AbelianGaugeCube& field, int direction, int layer, int time)
	{
		for (int n = 0; n < field.linear_size * field.linear_size; n++)
		{
			int n_shifted = n + cubic_size * time + planar_size * layer;
			int field_index = n_shifted * 4 + direction;
			assert(field_index < field.numVariables());
			const double theta = field.variable(field_index);

			site_pixels[n].SetFillColor(NormalMapYellow(theta, 0.3, 0.6, 0.8));
			//site_pixels[n].SetFillColor(GreenRedUniform(theta));
		}
	}
	void UpdateFluxColors(const AbelianGaugeCube& field, int layer, int time)
	{
		const std::vector<double> fluxes = field.getFluxes_z();
		for (int n = 0; n < planar_size; n++)
		{
			int n_shifted = n + cubic_size * time + planar_size * layer;
			site_pixels[n].SetFillColor(NormalMapYellow(fluxes[n_shifted], 0.3, 0.6, 0.8));
		}
	}
	sf::Color GreenRedUniform(const double& theta)
	{
		const int theta_c = int(std::abs(theta) * 510.0) % 510;
		sf::Uint8 r = sf::Uint8(255 - std::abs(theta_c - 255));
		sf::Uint8 g = sf::Uint8(std::abs(theta_c - 255));
		return sf::Color(r, g, 122u);
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
	int cubic_size = -1;
	int planar_size = -1;
	sf::RenderWindow& window;
	std::vector<Square> site_pixels;
	Vec2D offset{ 0.0f, 0.0f };
};