#include "MinMod.hpp"
#include <cmath>
#include <algorithm>

MinMod::MinMod(Boundary const& boundary):boundary_(boundary)
{}


MinMod::~MinMod()
{}

void MinMod::GetInterpolatedValues(vector<Primitive> const & cells, vector<double> const & edges, 
	vector<pair<Primitive, Primitive> >& values) const
{
	size_t N = edges.size();
	values.resize(N);
	Primitive slope;
	// Do bulk edges
	for (size_t i = 1; i < N - 2; ++i)
	{
		const Primitive sl = (cells[i] - cells[i - 1]) / (0.5*(edges[i + 1] - edges[i-1]));
		const Primitive sr = (cells[i+1] - cells[i]) / (0.5*(edges[i + 2] - edges[i]));
		const Primitive sc = (cells[i+1] - cells[i - 1]) / (0.5*(edges[i + 2]+edges[i + 1] - edges[i - 1]
			- edges[i]));
		if (sl.density*sr.density < 0)
			slope.density = 0;
		else
			slope.density = std::min(std::fabs(sl.density), std::min(std::fabs(sr.density), std::fabs(sc.density))) * (sl.density > 0 ?
				1 : -1);
		if (sl.pressure*sr.pressure < 0)
			slope.pressure = 0;
		else
			slope.pressure = std::min(std::fabs(sl.pressure), std::min(std::fabs(sr.pressure), std::fabs(sc.pressure))) * (sl.pressure > 0 ?
				1 : -1);
		if (sl.velocity*sr.velocity < 0)
			slope.velocity = 0;
		else
			slope.velocity = std::min(std::fabs(sl.velocity), std::min(std::fabs(sr.velocity), std::fabs(sc.velocity))) * (sl.velocity > 0 ?
				1 : -1);
		values[i].second = cells[i] - slope * (0.5*(edges[i + 1] - edges[i]));
		values[i+1].first = cells[i] + slope * (0.5*(edges[i + 1] - edges[i]));
	}
	// Do boundaries
	vector<Primitive> left = boundary_.GetBoundaryValues(cells, edges, 0);
	vector<Primitive> right = boundary_.GetBoundaryValues(cells, edges, edges.size()-1);
	values[0].first = left[0];
	values[0].second = left[1];
	values[1].first = left[2]; 
	values[N-2].second = right[0];
	values[N-1].first = right[1];
	values[N-1].second = right[2];
}
