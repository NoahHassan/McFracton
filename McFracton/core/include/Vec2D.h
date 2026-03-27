#pragma once

#include <cmath>

template<typename T>
class _Vec2D
{
public:
	_Vec2D() = default;
	constexpr _Vec2D(T x, T y) : x(x), y(y) {}
public:
	_Vec2D operator +(const _Vec2D& rhs) const
	{
		return { x + rhs.x, y + rhs.y };
	}
	_Vec2D& operator +=(const _Vec2D& rhs)
	{
		x += rhs.x;
		y += rhs.y;
		return *this;
	}
	_Vec2D operator -(const _Vec2D& rhs) const
	{
		return *this + (-rhs);
	}
	_Vec2D& operator -=(const _Vec2D& rhs)
	{
		return *this += (-rhs);
	}
	_Vec2D operator *(const T& rhs) const
	{
		return { x * rhs, y * rhs };
	}
	_Vec2D operator *=(const T& rhs)
	{
		x *= rhs;
		y *= rhs;
		return *this;
	}
	_Vec2D operator /(const T& rhs) const
	{
		return { x / rhs, y / rhs };
	}
	_Vec2D operator /=(const T& rhs)
	{
		x /= rhs;
		y /= rhs;
		return *this;
	}
	_Vec2D operator -() const
	{
		return { -x, -y };
	}
	bool operator ==(const _Vec2D& rhs) const
	{
		return x == rhs.x && y == rhs.y;
	}

#ifdef GUI_ENABLED
	operator sf::Vector2<T>() const
	{
		return sf::Vector2<T>(T(x), T(y));
	}
#endif
public:
	T GetLengthSqr() const
	{
		return x * x + y * y;
	}
	T GetLength() const
	{
		return sqrt(GetLengthSqr());
	}
	void ReflectAt(const _Vec2D& n)
	{
		*this -= Dot(*this, n) * T(2) * n;
	}
	_Vec2D Normalized() const
	{
		return *this / GetLength();
	}
public:
	static _Vec2D Zero()
	{
		return { T(0), T(0) };
	}
	static T Dot(_Vec2D v, _Vec2D w)
	{
		return v.x * w.x + v.y * w.y;
	}
	static T Cross(_Vec2D v, _Vec2D w)
	{
		return v.x * w.y - v.y * w.x;
	}
public:
	T x = T(0);
	T y = T(0);
};

template<typename T>
_Vec2D<T> operator *(const T& lhs, const _Vec2D<T>& rhs)
{
	return rhs * lhs;
}

template<typename T>
_Vec2D<T> Interpolate(const _Vec2D<T> a, const _Vec2D<T> b, float c)
{
	return a + (b - a) * c;
}

template<typename T>
static _Vec2D<T> _Vec2D<T>::Zero();

template<typename T>
static T _Vec2D<T>::Dot(_Vec2D<T> v, _Vec2D<T> w);

template<typename T>
static T _Vec2D<T>::Cross(_Vec2D<T> v, _Vec2D<T> w);

typedef _Vec2D<float> Vec2D;