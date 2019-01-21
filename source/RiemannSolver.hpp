#ifndef RS_HPP
#define RS_HPP 1

#include "Primitive.hpp"

struct RSsolution
{
	double velocity;
	double pressure;
};


class RiemannSolver
{
public:
	virtual ~RiemannSolver() {}

	virtual RSsolution Solve(Primitive const& left, Primitive const& right)const = 0;
};


#endif //RS_HPP