#include "Primitive.hpp"

Primitive::Primitive():density(0),pressure(0),velocity(0),entropy(0),energy(0),sticker(0),LastCool(0)
{}

Primitive::Primitive(Primitive const & other)
{
	density = other.density;
	pressure = other.pressure;
	velocity = other.velocity;
	entropy = other.entropy;
	energy = other.energy;
	sticker = other.sticker;
	LastCool = other.LastCool;
}

Primitive::Primitive(double Density, double Pressure, double Velocity, double Entropy, double Energy, 
	unsigned char Sticker):density(Density),pressure(Pressure),velocity(Velocity),entropy(Entropy),energy(Energy),
	sticker(Sticker),LastCool(0)
{}

Primitive operator-(Primitive const&left, Primitive const& right)
{
	return Primitive(left.density- right.density, left.pressure- right.pressure, left.velocity- right.velocity,
		left.entropy- right.entropy, left.energy- right.energy, left.sticker);
}

Primitive operator+(Primitive const&left, Primitive const& right)
{
	return Primitive(left.density + right.density, left.pressure + right.pressure, left.velocity + right.velocity,
		left.entropy + right.entropy, left.energy + right.energy, left.sticker);
}

Primitive operator*(Primitive const&p, double s)
{
	return Primitive(p.density*s,p.pressure*s,p.velocity*s,p.entropy*s,p.energy*s,p.sticker);
}

Primitive operator*(double s, Primitive const&p)
{
	return p * s;
}

Primitive operator/(Primitive const&p, double s)
{
	return Primitive(p.density/s, p.pressure/s, p.velocity/s,p.entropy/s,p.energy/s,p.sticker);
}

