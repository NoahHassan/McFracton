#pragma once

#include <SFML/Graphics.hpp>
#include <assert.h>

#include "Spiderweb.h"
#include "Square.h"

#define PI 3.1415926535897932384

class SpiderCanvas {
public:
	SpiderCanvas(sf::RenderWindow& window)
		:
		window(window)
	{}
	SpiderCanvas(sf::RenderWindow& window, Vec2D offset)
		:
		window(window), offset(offset)
	{}
	void Initialize(const Spiderweb& field, float square_size)
	{
		site_pixels = std::vector<Square>(field.linear_size * field.linear_size);
		localEnergies = std::vector<double>(field.linear_size * field.linear_size * field.temporal_size);
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
	void Draw(const Spiderweb& field, int direction)
	{
		UpdateFieldColors(field, direction);
		for (const auto& sq : site_pixels)
		{
			window.draw(sq);
		}
	}
	void Draw(const Spiderweb& field, int direction, int layer)
	{
		UpdateFieldColors(field, direction, layer);
		for (const auto& sq : site_pixels)
		{
			window.draw(sq);
		}
	}
	void Draw(const Spiderweb& field, std::vector<std::pair<std::vector<int>, int>> vortices, int direction, int layer)
	{
		UpdateFieldColors(field, direction, layer);
		//ColorVortices(vortices, layer * field.ss_size, (layer + 1) * field.ss_size);
		for (const auto& sq : site_pixels)
		{
			window.draw(sq);
		}
	}
	void DrawEnergy(const Spiderweb& field, int layer, double maxEnergy)
	{
		field.getEnergy(localEnergies);
		UpdateEnergyColors(field, layer, maxEnergy);
		for (const auto& sq : site_pixels)
		{
			window.draw(sq);
		}
	}
private:
	void UpdateFieldColors(const Spiderweb& field, int direction)
	{
		for (int n = 0; n < field.linear_size * field.linear_size; n++)
		{
			int field_index = n * 3 + direction;
			assert(field_index < field.n_site_variables);
			const double& theta = field.getSite(field_index);

			//site_pixels[n].SetFillColor(NormalMapYellow(theta, 0.3, 0.6, 0.8));
			site_pixels[n].SetFillColor(GreenRedUniform(theta));
		}
	}
	void UpdateFieldColors(const Spiderweb& field, int direction, int layer)
	{
		for (int n = 0; n < field.linear_size * field.linear_size; n++)
		{
			int n_shifted = n + (int)site_pixels.size() * layer;
			int field_index = n_shifted * 3 + direction;
			assert(field_index < field.n_site_variables);
			const double& theta = field.getSite(field_index);

			site_pixels[n].SetFillColor(NormalMapYellow(theta, 0.3, 0.6, 0.8));
			//site_pixels[n].SetFillColor(GreenRedUniform(theta));
		}
	}
	void UpdateEnergyColors(const Spiderweb& field, int layer, double maxEnergy)
	{
		for (int n = 0; n < field.linear_size * field.linear_size; n++)
		{
			int n_shifted = n + (int)site_pixels.size() * layer;
			//site_pixels[n].SetFillColor(NormalMapYellow(localEnergies[n_shifted], 0.3, 0.6, 0.8));
			//site_pixels[n].SetFillColor(GreenRedUniform(localEnergies[n_shifted]));
			site_pixels[n].SetFillColor(BlackWhite(localEnergies[n_shifted] / maxEnergy));
		}
	}
	void ColorVortices(std::vector<std::pair<std::vector<int>, int>> vortices, int range_min, int range_max)
	{
		for (auto it = vortices.begin(); it != vortices.end(); ++it)
		{
			auto vortex = (*it).first;
			int type = (*it).second;

			for (int n = 0; n < vortex.size(); n++)
			{
				int vortex_index = vortex[n];
				int num = 0;
				if (vortex_index > range_min && vortex_index < range_max)
				{
					if (type == 1) {
						num++;
						site_pixels[vortex_index - range_min].SetFillColor(sf::Color::Blue);
					}
					else {
						site_pixels[vortex_index - range_min].SetFillColor(sf::Color::Yellow);
					}
				}
				if (num == 1) {
					const auto meh = 1 + 1;
				}
			}
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
		// map cos in range (-3,3) to range (0,2) and apply triangle function
		const double mapped_cos = abs(cosine / 2.0 + 0.5); // lowest energies colored black, highest white
		//const double n_val = -abs(mapped_cos - 1.0) + 1.0;
		const int color_val = int(mapped_cos * 510) % 510;
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
		if (n == 0)
			return sf::Color::White;
		if (n > 0)
			return sf::Color::Blue;
		if (n < 0)
			return sf::Color::Red;
	}
private:
	sf::RenderWindow& window;
	std::vector<Square> site_pixels;
	std::vector<double> localEnergies;
	Vec2D offset{ 0.0f, 0.0f };
};