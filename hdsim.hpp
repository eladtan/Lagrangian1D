#ifndef HDSIM_HPP
#define HDSIM_HPP 1

#include "Primitive.hpp"
#include "Interpolation.hpp"
#include "ideal_gas.hpp"
#include "ExactRS.hpp"
#include "Extensive.hpp"
#include "SourceTerm.hpp"
#include "Boundary.hpp"
#include "Geometry.hpp"
#include <vector>

using namespace std;

class hdsim
{
private:
	const double cfl_;
	vector<Primitive> cells_;
	vector<double> edges_;
	Interpolation const& interpolation_;
	IdealGas const& eos_;
	ExactRS const& rs_;
	double time_;
	size_t cycle_;
	vector<pair<Primitive, Primitive> > interp_values_;
	vector<RSsolution> rs_values_;
	vector<Extensive> extensives_;
	SourceTerm const& source_;
	Geometry const& geo_;
	BoundarySolution const* BoundarySolution_;
public:
	hdsim(double cfl,vector<Primitive> const& cells,vector<double> const& edges,Interpolation const& interp,
		IdealGas const& eos,ExactRS const& rs,SourceTerm const& source,Geometry const& geo, BoundarySolution const* BS=0);
	~hdsim();
	void TimeAdvance2();
	void TimeAdvance();
	double GetTime()const;
	vector<Primitive>const& GetCells()const;
	vector<Primitive>& GetCells();
	vector<Extensive>const& GetExtensives()const;
	vector<Extensive>& GetExtensives();
	vector<double> const& GetEdges()const;
	size_t GetCycle()const;
	void SetTime(double t);
	void SetCycle(size_t cyc);
	void ReCalcCells(vector<Extensive> const& extensives);
};
#endif //HDSIM_HPP
