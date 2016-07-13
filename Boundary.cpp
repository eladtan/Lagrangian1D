#include "Boundary.hpp"
#include <algorithm>
#include<cmath>

BoundarySolution::~BoundarySolution() {}

VacuumInFlow::VacuumInFlow(bool calc_left, bool calc_right) :calc_left_(calc_left), calc_right_(calc_right_) {}

pair<RSsolution, RSsolution> VacuumInFlow::GetBoundaryValues(vector<Primitive> const& cells)const
{
	RSsolution left, right;
	left.pressure = 0;
	right.pressure = 0;
	left.velocity = cells[0].velocity;
	right.velocity = cells.back().velocity;
	return pair<RSsolution, RSsolution>(left, right);
}

std::pair<bool, bool> VacuumInFlow::ShouldCalc()const
{
	return pair<bool, bool>(calc_left_, calc_right_);
}

Boundary::~Boundary()
{}

vector<Primitive> RigidWall::GetBoundaryValues(vector<Primitive> const & cells,
	vector<double> const & edges, size_t index,double time) const
{
	vector<Primitive> res(3);
	if (index == 0)
	{
		Primitive slope = (cells[1] - cells[0]) / (0.5*(edges[2] - edges[0]));
		slope.density = 0;
		slope.entropy = 0;
		//slope.velocity = 0;
		slope.pressure = 0;
		res[1] = cells[0] - slope*(0.5*(edges[1] - edges[0]));
		res[0] = res[1];
		res[0].velocity = -res[1].velocity;
		res[2] = cells[0]+slope*(0.5*(edges[1] - edges[0]));
	}
	else
	{
		size_t N = edges.size();
		Primitive slope = (cells[N-2] - cells[N-3]) / (0.5*(edges[N-1] - edges[N-3]));
		slope.density = 0;
		slope.entropy = 0;
		//slope.velocity = 0;
		slope.pressure = 0;
		res[1] = cells[N-2] + slope*(0.5*(edges[N-1] - edges[N-2]));
		res[2] = res[1];
		res[2].velocity = -res[1].velocity;
		res[0] = cells[N-2] - slope*(0.5*(edges[N-1] - edges[N-2]));
	}
	return res;
}

vector<Primitive> FreeFlow::GetBoundaryValues(vector<Primitive> const & cells, vector<double> const & edges,
	size_t index, double time) const
{
	vector<Primitive> res(3);
	if (index == 0)
	{
		Primitive slope = (cells[1] - cells[0]) / (0.5*(edges[2] - edges[0]));
		res[1] = cells[0] - slope*(0.5*(edges[1] - edges[0]));
		res[0] = res[1];
		res[2] = cells[0] + slope*(0.5*(edges[1] - edges[0]));
	}
	else
	{
		size_t N = edges.size();
		Primitive slope = (cells[N - 2] - cells[N - 3]) / (0.5*(edges[N - 1] - edges[N - 3]));
		res[1] = cells[N - 2] + slope*(0.5*(edges[N - 1] - edges[N - 2]));
		res[2] = res[1];
		res[0] = cells[N - 2] - slope*(0.5*(edges[N - 1] - edges[N - 2]));
	}
	return res;
}

vector<Primitive> Periodic::GetBoundaryValues(vector<Primitive> const & cells, vector<double> const & edges, size_t index,
	double time) const
{
	vector<Primitive> res(3);
	size_t N = edges.size();
	double L = edges[N - 1] - edges[0];
	Primitive slope0,slopeN;
	Primitive sr = (cells[1] - cells[0]) / (0.5*(edges[2] - edges[0]));
	Primitive sl = (cells[0] - cells[N - 2]) / (0.5*(edges[1] - edges[N - 2] + L));
	Primitive sc = (cells[1] - cells[N - 2]) / (0.5*(edges[2] + edges[1] - edges[0] - edges[N - 2] + L));
	if (sl.density*sr.density < 0)
		slope0.density = 0;
	else
		slope0.density = min(fabs(sl.density), min(fabs(sr.density), fabs(sc.density))) * (sl.density > 0 ?
			1 : -1);
	if (sl.pressure*sr.pressure < 0)
		slope0.pressure = 0;
	else
		slope0.pressure = min(fabs(sl.pressure), min(fabs(sr.pressure), fabs(sc.pressure))) * (sl.pressure > 0 ?
			1 : -1);
	if (sl.velocity*sr.velocity < 0)
		slope0.velocity = 0;
	else
		slope0.velocity = min(fabs(sl.velocity), min(fabs(sr.velocity), fabs(sc.velocity))) * (sl.velocity > 0 ?
			1 : -1);
	sr = sl;
	sl = (cells[N - 2] - cells[N - 3]) / (0.5*(edges[N - 1] - edges[N - 3]));
	sc = (cells[0] - cells[N - 3]) / (0.5*(edges[N - 1] + edges[0] - edges[N - 3] - edges[N - 2] + L));
	if (sl.density*sr.density < 0)
		slopeN.density = 0;
	else
		slopeN.density = min(fabs(sl.density), min(fabs(sr.density), fabs(sc.density))) * (sl.density > 0 ?
			1 : -1);
	if (sl.pressure*sr.pressure < 0)
		slopeN.pressure = 0;
	else
		slopeN.pressure = min(fabs(sl.pressure), min(fabs(sr.pressure), fabs(sc.pressure))) * (sl.pressure > 0 ?
			1 : -1);
	if (sl.velocity*sr.velocity < 0)
		slopeN.velocity = 0;
	else
		slopeN.velocity = min(fabs(sl.velocity), min(fabs(sr.velocity), fabs(sc.velocity))) * (sl.velocity > 0 ?
			1 : -1);

	if (index == 0)
	{
		res[1] = cells[0] - slope0*(0.5*(edges[1] - edges[0]));
		res[0] = cells[N-2]+slopeN*(0.5*(edges[N-1] - edges[N-2]));
		res[2] = cells[0] + slope0*(0.5*(edges[1] - edges[0]));
	}
	else
	{
		res[1] = cells[N - 2] + slopeN*(0.5*(edges[N - 1] - edges[N - 2]));
		res[2] = cells[0] - slope0*(0.5*(edges[1] - edges[0]));
		res[0] = cells[N - 2] - slopeN*(0.5*(edges[N - 1] - edges[N - 2]));
	}
	return res;
}

SeveralBoundary::SeveralBoundary(Boundary const & left, Boundary const & right):left_(left),right_(right)
{}

vector<Primitive> SeveralBoundary::GetBoundaryValues(vector<Primitive> const & cells, vector<double> const & edges, size_t index, 
	double time) const
{
	if(index==0)
		return left_.GetBoundaryValues(cells, edges, index,time);
	else
		return right_.GetBoundaryValues(cells, edges, index,time);
}

ConstantPrimitive::ConstantPrimitive(Primitive outer):outer_(outer)
{}

vector<Primitive> ConstantPrimitive::GetBoundaryValues(vector<Primitive> const & cells, vector<double> const & edges, size_t index, 
	double time) const
{
	vector<Primitive> res(3);
	size_t N = edges.size();
	Primitive left, center, right;
	if (index == 0)
	{
		left = outer_;
		center = cells[0];
		right = cells[1];
	}
	else
	{
		left = cells[N-3];
		center = cells[N-2];
		right = outer_;
	}
	double dxl, dxr, dxc;
	if (index == 0)
	{
		dxl = edges[1] - edges[0];
		dxr = 0.5*(edges[2] - edges[0]);
		dxc = 0.5*edges[2] + edges[1] - 1.5*edges[0];
	}
	else
	{
		dxl = 0.5*(edges[N-1] - edges[N-3]);
		dxr = edges[N-1] - edges[N-2];
		dxc = 1.5*edges[N - 1] - edges[N - 2] - 0.5*edges[N - 3];
	}
	Primitive sr = (right - center) / dxr;
	Primitive sl = (center - left) / dxl;
	Primitive sc = (right - left) / dxc;
	Primitive slope0;
	if (sl.density*sr.density < 0)
		slope0.density = 0;
	else
		slope0.density = min(fabs(sl.density), min(fabs(sr.density), fabs(sc.density))) * (sl.density > 0 ?
			1 : -1);
	if (sl.pressure*sr.pressure < 0)
		slope0.pressure = 0;
	else
		slope0.pressure = min(fabs(sl.pressure), min(fabs(sr.pressure), fabs(sc.pressure))) * (sl.pressure > 0 ?
			1 : -1);
	if (sl.velocity*sr.velocity < 0)
		slope0.velocity = 0;
	else
		slope0.velocity = min(fabs(sl.velocity), min(fabs(sr.velocity), fabs(sc.velocity))) * (sl.velocity > 0 ?
			1 : -1);
	if (index == 0)
	{
		res[1] = cells[0] - slope0*(0.5*(edges[1] - edges[0]));
		res[0] = outer_;
		res[2] = cells[0] + slope0*(0.5*(edges[1] - edges[0]));
	}
	else
	{
		res[1] = cells[N - 2] + slope0*(0.5*(edges[N - 1] - edges[N - 2]));
		res[2] = outer_;
		res[0] = cells[N - 2] - slope0*(0.5*(edges[N - 1] - edges[N - 2]));
	}
	return res;
}
