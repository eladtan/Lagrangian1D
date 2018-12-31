#ifndef MPI_COMM_HPP
#define MPI_COMM_HPP 1
#ifdef RICH_MPI
#include <mpi.h>
#include "Extensive.hpp"
#include "Primitive.hpp"
#include <vector>

std::array<Primitive,4> SendRecvPrimitive(std::vector<Primitive> const& cells);

std::array<double, 4> SendRecvEdges(std::vector<double> const& edges);

void RedistributeExtensives(std::vector<Extensive> &cells,std::vector<double> &edges);

void ConsolidateData(std::vector<Primitive> &cells, std::vector<double> &edges,std::vector<std::vector<
	double> > & append);

#endif

#endif //MPI_COMM_HPP