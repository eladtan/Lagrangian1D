#ifndef HLLC_HPP
#define HLLC_HPP 1

#include "RiemannSolver.hpp"
#include "ideal_gas.hpp"

class Hllc :public RiemannSolver
{
private:
	IdealGas const& eos_;
	const bool iter_;

public:
	Hllc(IdealGas const& eos,bool iter);
	~Hllc();

	RSsolution Solve(Primitive const& left, Primitive const& right)const;
};


#endif //HLLC_HPP