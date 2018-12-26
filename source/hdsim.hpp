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
	double TotalEcool_;
	vector<pair<Primitive, Primitive> > interp_values_;
	vector<RSsolution> rs_values_;
	vector<Extensive> extensives_;
	SourceTerm const& source_;
	Geometry const& geo_;
	const double AMR_ratio_;
	BoundarySolution const* BoundarySolution_;
	double dt_suggest_;
	void AMR(void);
public:
	hdsim(double cfl,vector<Primitive> const& cells,vector<double> const& edges,Interpolation const& interp,
		IdealGas const& eos,ExactRS const& rs,SourceTerm const& source,Geometry const& geo, const double AMR_ratio = 0, 
		BoundarySolution const* BS=0);
	~hdsim();
	void TimeAdvance2();
	void TimeAdvance();
	double GetTime()const;
	double GetEcool()const;
	vector<Primitive>const& GetCells()const;
	vector<Primitive>& GetCells();
	vector<Extensive>const& GetExtensives()const;
	vector<Extensive>& GetExtensives();
	vector<double> const& GetEdges()const;
	vector<double>& GetEdges();
	size_t GetCycle()const;
	void SetTime(double t);
	void SetEcool(double E);
	void SetCycle(size_t cyc);
	void ReCalcCells(vector<Extensive> & extensives);
	void ReCalcExtensives(vector<Primitive> const& cells);
	void SuggestTimeStep(double dt);
};
#endif //HDSIM_HPP
