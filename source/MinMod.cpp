#include "MinMod.hpp"
#include <cmath>
#include <algorithm>
#ifdef RICH_MPI
#include <array>
#include "mpi_comm.hpp"
#endif

namespace
{
	Primitive GetSlope(Primitive const& left, Primitive const& center, Primitive const& right,
		double e0, double e1, double e2, double e3)
	{
		Primitive slope;
		const Primitive sl = (center - left) / (0.5*(e2 - e0));
		const Primitive sr = (right - center) / (0.5*(e3 - e1));
		const Primitive sc = (right - left) / (0.5*(e3 + e2 - e0 - e1));
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
		return slope;
	}
}

MinMod::MinMod(Boundary const& boundary):boundary_(boundary)
{}


MinMod::~MinMod()
{}

void MinMod::GetInterpolatedValues(vector<Primitive> const & cells, vector<double> const & edges, 
	vector<pair<Primitive, Primitive> >& values,double time) const
{
	size_t N = edges.size();
	values.resize(N);
	Primitive slope;
	// Do bulk edges
	for (size_t i = 1; i < N - 2; ++i)
	{
		slope = GetSlope(cells[i - 1], cells[i], cells[i + 1], edges[i - 1], edges[i], 
			edges[i + 1], edges[i + 2]);
		values[i].second = cells[i] - slope * (0.5*(edges[i + 1] - edges[i]));
		values[i+1].first = cells[i] + slope * (0.5*(edges[i + 1] - edges[i]));
	}
	// Do boundaries
#ifdef RICH_MPI
	std::array<Primitive, 4> ghost_cells = SendRecvPrimitive(cells);
	std::array<double, 4> ghost_edges = SendRecvEdges(edges);
	int rank = 0, ws = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	if (rank == 0)
	{
		slope = GetSlope(cells[N - 3], cells[N - 2], ghost_cells[2], edges[N - 3], edges[N - 2],
			edges[N - 1], ghost_edges[2]);
		values[N - 1].first = cells[N - 2] + slope * 0.5*(edges[N - 1] - edges[N - 2]);
		values[N - 2].second = cells[N - 2] - slope * 0.5*(edges[N - 1] - edges[N - 2]);
		slope = GetSlope(cells[N - 2], ghost_cells[2], ghost_cells[3], edges[N - 2], edges[N - 1],
			ghost_edges[2], ghost_edges[3]);
		values[N - 1].second = ghost_cells[2] - slope * 0.5 * (ghost_edges[2] - edges[N - 1]);
		vector<Primitive> left = boundary_.GetBoundaryValues(cells, edges, 0, time);
		values[0].first = left[0];
		values[0].second = left[1];
		values[1].first = left[2];
	}
	else
	{
		if (rank == (ws - 1))
		{
			slope = GetSlope(ghost_cells[1], cells[0], cells[1], ghost_edges[1], edges[0],
				edges[1], edges[2]);
			values[1].first = cells[0] + slope * 0.5*(edges[1] - edges[0]);
			values[0].second = cells[0] - slope * 0.5*(edges[1] - edges[0]);
			slope = GetSlope(ghost_cells[0], ghost_cells[1], cells[0], ghost_edges[0], ghost_edges[1],
				edges[0], edges[1]);
			values[0].first = ghost_cells[1] + slope * 0.5 * (edges[0] - ghost_edges[1]);
			vector<Primitive> right = boundary_.GetBoundaryValues(cells, edges, edges.size() - 1, time);
			values[N - 2].second = right[0];
			values[N - 1].first = right[1];
			values[N - 1].second = right[2];
		}
		else
		{
			slope = GetSlope(ghost_cells[1], cells[0], cells[1], ghost_edges[1], edges[0],
				edges[1], edges[2]);
			values[1].first = cells[0] + slope * 0.5*(edges[1] - edges[0]);
			values[0].second = cells[0] - slope * 0.5*(edges[1] - edges[0]);
			slope = GetSlope(ghost_cells[0], ghost_cells[1], cells[0], ghost_edges[0], ghost_edges[1],
				edges[0], edges[1]);
			values[0].first = ghost_cells[1] + slope * 0.5 * (edges[0] - ghost_edges[1]);

			slope = GetSlope(cells[N-3], cells[N-2], ghost_cells[2], edges[N-3],edges[N-2],
				edges[N-1], ghost_edges[2]);
			values[N-1].first = cells[N-2] + slope * 0.5*(edges[N-1] - edges[N-2]);
			values[N-2].second = cells[N-2] - slope * 0.5*(edges[N-1] - edges[N-2]);
			slope = GetSlope(cells[N - 2], ghost_cells[2], ghost_cells[3], edges[N - 2], edges[N - 1],
				ghost_edges[2], ghost_edges[3]);
			values[N-1].second = ghost_cells[2] - slope * 0.5 * (ghost_edges[2] - edges[N-1]);
		}
	}
#else
	vector<Primitive> left = boundary_.GetBoundaryValues(cells, edges, 0, time);
	values[0].first = left[0];
	values[0].second = left[1];
	values[1].first = left[2];
	vector<Primitive> right = boundary_.GetBoundaryValues(cells, edges, edges.size() - 1, time);
	values[N - 2].second = right[0];
	values[N - 1].first = right[1];
	values[N - 1].second = right[2];
#endif
}
