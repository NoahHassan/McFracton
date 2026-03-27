#pragma once

#include <memory>

#include "Shape.h"
#include "Vec2D.h"

class Square : public Shape
{
public:
	Square() = default;
	Square(Vec2D center, float size)
	{
		sf::RectangleShape rect(Vec2D{ size, size });
		rect.setOrigin(Vec2D{ size / 2.0f, size / 2.0f });
		shape = std::make_shared<sf::RectangleShape>(rect);
		SetPosition(center);
	}
	virtual ~Square() = default;
};