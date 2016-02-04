#ifndef HDSIM_HPP
#define HDSIM_HPP 1

#include "Primitive.hpp"
#include "MinMod.hpp"
#include "ideal_gas.hpp"
#include "ExactRS.hpp"
#include "Extensive.hpp"
#include "SourceTerm.hpp"
#include <vector>

using namespace std;

class hdsim
{
private:
	const double cfl_;
	vector<Primitive> cells_;
	vector<double> edges_;
	MinMod const& interpolation_;
	IdealGas const& eos_;
	ExactRS const& rs_;
	double time_;
	size_t cycle_;
	vector<pair<Primitive, Primitive> > interp_values_;
	vector<RSsolution> rs_values_;
	vector<Extensive> extensives_;
	SourceTerm const& source_;
public:
	hdsim(double cfl,vector<Primitive> const& cells,vector<double> const& edges,MinMod const& interp,
		IdealGas const& eos,ExactRS const& rs,SourceTerm const& source);
	~hdsim();
	void TimeAdvance2();
	double GetTime()const;
	vector<Primitive>const& GetCells()const;
	vector<double> const& GetEdges()const;
	size_t GetCycle()const;
	void SetTime(double t);
};
#endif //HDSIM_HPP
