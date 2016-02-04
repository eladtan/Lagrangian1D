#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP 1

#include "Primitive.hpp"
#include <vector>

using namespace std;

class Boundary
{
public:
	virtual ~Boundary();

	virtual vector<Primitive> GetBoundaryValues(vector<Primitive> const& cells, vector<double> const&
		edges,size_t index)const=0;
};

class RigidWall : public Boundary
{
public:
	vector<Primitive> GetBoundaryValues(vector<Primitive> const& cells, vector<double> const&
		edges, size_t index)const;
};

class FreeFlow : public Boundary
{
public:
	vector<Primitive> GetBoundaryValues(vector<Primitive> const& cells, vector<double> const&
		edges, size_t index)const;
};

class Periodic : public Boundary
{
public:
	vector<Primitive> GetBoundaryValues(vector<Primitive> const& cells, vector<double> const&
		edges, size_t index)const;
};
#endif //BOUNDARY_HPP
