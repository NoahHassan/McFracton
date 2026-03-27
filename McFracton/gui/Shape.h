#pragma once

#define GUI_ENABLED

#include <SFML/Graphics.hpp>
#include <memory>

#include "Vec2D.h"

class Shape : public sf::Drawable
{
public:
	virtual ~Shape() {};
public:
	virtual void draw(sf::RenderTarget& target, sf::RenderStates states) const;
public:
	void SetPosition(Vec2D position);
	void SetScale(Vec2D scale);
	void SetRotation(float rotation);
	void SetFillColor(sf::Color color);
	void SetOutlineColor(sf::Color color);
	void SetOutlineThickness(float thickness);
protected:
	std::shared_ptr<sf::Shape> shape;
private:
	Vec2D position;
	Vec2D scale;
	float rotation = 0.0f;
};