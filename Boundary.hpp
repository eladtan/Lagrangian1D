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

class ConstantPrimitive : public Boundary
{
private:
	Primitive outer_;
public:
	ConstantPrimitive(Primitive outer);

	vector<Primitive> GetBoundaryValues(vector<Primitive> const& cells, vector<double> const&
		edges, size_t index)const;
};

class SeveralBoundary : public Boundary
{
private:
	Boundary const& left_;
	Boundary const& right_;
public:
	SeveralBoundary(Boundary const& left, Boundary const& right);

	vector<Primitive> GetBoundaryValues(vector<Primitive> const& cells, vector<double> const&
		edges, size_t index)const;

};
#endif //BOUNDARY_HPP
