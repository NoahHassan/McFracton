#pragma once

#include <SFML/Graphics.hpp>

#include "XYSquare.h"
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
	void Initialize(const XYSquare& field, float square_size)
	{
		site_pixels = std::vector<Square>(field.nSites);
		for (int n = 0; n < field.nSites; n++)
		{
			int ny = n / field.length;
			int nx = n - ny * field.length;

			Square sq(
				Vec2D(
					square_size * (nx - field.length / 2) + window.getSize().x / 2, 
					square_size * (ny - field.length / 2) + window.getSize().y / 2
				) + offset,
				square_size
			);

			sq.SetOutlineThickness(2.0);
			sq.SetOutlineColor(sf::Color::Black);
			
			site_pixels[n] = std::move(sq);
		}
	}
	void Draw(const XYSquare& field)
	{
		UpdateFieldColors(field);
		for (const auto& sq : site_pixels)
		{
			window.draw(sq);
		}
	}
	void Draw(const XYSquare& field, std::vector<std::pair<std::vector<int>, int>> vortices)
	{
		UpdateFieldColors(field);
		ColorVortices(vortices);
		for (const auto& sq : site_pixels)
		{
			window.draw(sq);
		}
	}
private:
	void UpdateFieldColors(const XYSquare& field)
	{
		// There are two points that are equivalent here
		// so may need to change this color wheel to something better
		for (int n = 0; n < field.nSites; n++)
		{
			const double& theta = field.getSite(n);

			site_pixels[n].SetFillColor(NormalMapYellow(theta, 0.3, 0.6, 0.8));
			//site_pixels[n].SetFillColor(GreenRedUniform(theta));
		}
	}
	void ColorVortices(std::vector<std::pair<std::vector<int>, int>> vortices)
	{
		std::for_each(vortices.begin(), vortices.end(),
			[&](std::pair<std::vector<int>, int> v) {
				std::for_each(v.first.begin(), v.first.end(),
				[&](int i) {
						if (v.second == 1)
							site_pixels[i].SetFillColor(sf::Color::Blue);
						else if (v.second == -1)
							site_pixels[i].SetFillColor(sf::Color::Yellow);
					}
				);
			}
		);
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
private:
	sf::RenderWindow& window;
	std::vector<Square> site_pixels;
	Vec2D offset{ 0.0f, 0.0f };
};