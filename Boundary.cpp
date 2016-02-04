#include "Boundary.hpp"
#include <algorithm>

Boundary::~Boundary()
{}

vector<Primitive> RigidWall::GetBoundaryValues(vector<Primitive> const & cells,
	vector<double> const & edges, size_t index) const
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
	size_t index) const
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

vector<Primitive> Periodic::GetBoundaryValues(vector<Primitive> const & cells, vector<double> const & edges, size_t index) const
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
		slope0.density = min(abs(sl.density), min(abs(sr.density), abs(sc.density))) * (sl.density > 0 ?
			1 : -1);
	if (sl.pressure*sr.pressure < 0)
		slope0.pressure = 0;
	else
		slope0.pressure = min(abs(sl.pressure), min(abs(sr.pressure), abs(sc.pressure))) * (sl.pressure > 0 ?
			1 : -1);
	if (sl.velocity*sr.velocity < 0)
		slope0.velocity = 0;
	else
		slope0.velocity = min(abs(sl.velocity), min(abs(sr.velocity), abs(sc.velocity))) * (sl.velocity > 0 ?
			1 : -1);
	sr = sl;
	sl = (cells[N - 2] - cells[N - 3]) / (0.5*(edges[N - 1] - edges[N - 3]));
	sc = (cells[0] - cells[N - 3]) / (0.5*(edges[N - 1] + edges[0] - edges[N - 3] - edges[N - 2] + L));
	if (sl.density*sr.density < 0)
		slopeN.density = 0;
	else
		slopeN.density = min(abs(sl.density), min(abs(sr.density), abs(sc.density))) * (sl.density > 0 ?
			1 : -1);
	if (sl.pressure*sr.pressure < 0)
		slopeN.pressure = 0;
	else
		slopeN.pressure = min(abs(sl.pressure), min(abs(sr.pressure), abs(sc.pressure))) * (sl.pressure > 0 ?
			1 : -1);
	if (sl.velocity*sr.velocity < 0)
		slopeN.velocity = 0;
	else
		slopeN.velocity = min(abs(sl.velocity), min(abs(sr.velocity), abs(sc.velocity))) * (sl.velocity > 0 ?
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
