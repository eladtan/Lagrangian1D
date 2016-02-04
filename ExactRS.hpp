#ifndef EXACTRS_HPP
#define EXACTRS_HPP 1

#include "Primitive.hpp"

struct RSsolution
{
	double velocity;
	double pressure;
};


class ExactRS
{
private:
	const double gamma_;

public:
	ExactRS(double gama);
	~ExactRS();

	RSsolution Solve(Primitive const& left, Primitive const& right)const;
};


#endif //EXACTRS_HPP