#include "Primitive.hpp"

Primitive::Primitive():density(0),pressure(0),velocity(0),entropy(0),energy(0)
{}

Primitive::Primitive(Primitive const & other)
{
	density = other.density;
	pressure = other.pressure;
	velocity = other.velocity;
	entropy = other.entropy;
	energy = other.energy;
}

Primitive::Primitive(double Density, double Pressure, double Velocity, double Entropy, double Energy):density(Density),
	pressure(Pressure),velocity(Velocity),entropy(Entropy),energy(Energy)
{}

Primitive Primitive::operator-(Primitive const & other)const
{
	return Primitive(this->density-other.density,this->pressure-other.pressure,this->velocity-other.velocity,
		this->entropy-other.entropy,this->energy-other.energy);
}

Primitive Primitive::operator+(Primitive const & other)const
{
	return Primitive(this->density + other.density, this->pressure + other.pressure, this->velocity + other.velocity,
		this->entropy+other.entropy,this->energy+other.energy);
}

Primitive Primitive::operator*(double s)
{
	return Primitive(this->density*s,this->pressure*s,this->velocity*s,this->entropy*s,this->energy*s);
}

Primitive Primitive::operator/(double s)
{
	return Primitive(this->density/s, this->pressure/s, this->velocity/s,this->entropy/s,this->energy/s);
}

