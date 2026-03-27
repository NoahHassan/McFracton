#include "Shape.h"

void Shape::draw(sf::RenderTarget& target, sf::RenderStates states) const
{
	target.draw(*shape, states);
}

void Shape::SetPosition(Vec2D position)
{
	shape->setPosition(position);
}

void Shape::SetScale(Vec2D scale)
{
	shape->setScale(scale);
}

void Shape::SetRotation(float rotation)
{
	shape->setRotation(rotation);
}

void Shape::SetFillColor(sf::Color color)
{
	shape->setFillColor(color);
}

void Shape::SetOutlineColor(sf::Color color)
{
	shape->setOutlineColor(color);
}

void Shape::SetOutlineThickness(float thickness)
{
	shape->setOutlineThickness(thickness);
}
