#ifndef PRIMITIVE_HPP
#define PRIMITIVE_HPP 1

class Primitive
{
public:
	double density;
	double pressure;
	double velocity;
	double entropy;
	double energy;
	unsigned char sticker;
	double LastCool;

	Primitive();

	Primitive(Primitive const& other);

	Primitive(double Density, double Pressure, double Velocity,double Entropy,double Energy,unsigned char Sticker);
};


Primitive operator-(Primitive const&left, Primitive const& right);

Primitive operator+(Primitive const&left, Primitive const& right);

Primitive operator*(Primitive const&p, double s);

Primitive operator*(double s,Primitive const&p);

Primitive operator/(Primitive const&p,  double s);

#endif //PRIMITIVE_HPP
