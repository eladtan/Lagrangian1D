#include "hdsim.hpp"
#include <algorithm>


hdsim::hdsim(double cfl, vector<Primitive> const& cells, vector<double> const& edges, MinMod const& interp,
	IdealGas const& eos, ExactRS const& rs,SourceTerm const& source):cfl_(cfl),cells_(cells),edges_(edges),interpolation_(interp),eos_(eos),
	rs_(rs),time_(0),cycle_(0),extensives_(vector<Extensive>()),source_(source)
{
	size_t N = cells.size();
	extensives_.resize(N);
	for (size_t i = 0; i < N; ++i)
	{
		double vol = edges[i + 1] - edges[i];
		extensives_[i].mass = cells[i].density*vol;
		extensives_[i].momentum = extensives_[i].mass*cells[i].velocity;
		extensives_[i].energy = 0.5*extensives_[i].momentum*extensives_[i].momentum / extensives_[i].mass +
			eos_.dp2e(cells[i].density, cells[i].pressure)*extensives_[i].mass;
	}
}


hdsim::~hdsim()
{}

namespace
{
	double GetTimeStep(vector<Primitive> const& cells, vector<double> const& edges, IdealGas const& eos, double cfl)
	{
		double dt = (edges[1] - edges[0]) / eos.dp2c(cells[0].density, cells[0].pressure);
		size_t N = cells.size();
		for (size_t i = 1; i < N;++i)
			dt = min(dt,(edges[i+1] - edges[i]) / eos.dp2c(cells[i].density, cells[i].pressure));
		return dt*cfl;
	}

	void GetRSvalues(vector<pair<Primitive,Primitive> > const& interp_values, ExactRS const&rs,
		vector<RSsolution> &res)
	{
		size_t N = interp_values.size();
		res.resize(N);
		for (size_t i = 0; i < N; ++i)
			res[i] = rs.Solve(interp_values[i].first, interp_values[i].second);
	}

	void UpdateExtensives(vector<Extensive> &cells, vector<RSsolution> const& rs_values_,double dt)
	{
		size_t N = cells.size();
		for (size_t i = 0; i < N; ++i)
		{
			cells[i].momentum -= (rs_values_[i + 1].pressure - rs_values_[i].pressure)*dt;
			cells[i].energy -= (rs_values_[i + 1].pressure*rs_values_[i + 1].velocity - rs_values_[i].pressure*
				rs_values_[i].velocity)*dt;
		}
	}

	void UpdateEdges(vector<double> &edges, vector<RSsolution> const& rs_values_,double dt)
	{
		size_t N = edges.size();
		for (size_t i = 0; i < N; ++i)
			edges[i] += rs_values_[i].velocity*dt;
	}

	bool ShouldUseEntropy(Primitive const& cell, vector<RSsolution> const& rsvalues, size_t index)
	{
		double ek = cell.velocity*cell.velocity;
		double et = cell.pressure / cell.density;
		if (et > 0.01*ek)
			return false;
		double dv = rsvalues[index + 1].velocity - rsvalues[index].velocity;
		if (dv > 0)
			return true;
		if (dv*dv > 0.1*et)
			return false;
		else
			return true;
	}

	void UpdateCells(vector<Extensive> &extensive, vector<double> const& edges, IdealGas const& eos,
		vector<Primitive> &cells,vector<RSsolution> const& rsvalues)
	{
		size_t N = cells.size();
		for (size_t i = 0; i < N; ++i)
		{
			double vol = edges[i + 1] - edges[i];
			cells[i].density = extensive[i].mass / vol;
			cells[i].velocity = extensive[i].momentum / extensive[i].mass;
			if (ShouldUseEntropy(cells[i], rsvalues, i))
				cells[i].pressure = eos.sd2p(cells[i].entropy, cells[i].density);
			else
				cells[i].pressure = eos.de2p(cells[i].density, (extensive[i].energy - 0.5*extensive[i].momentum*extensive[i].momentum
					/ extensive[i].mass) / extensive[i].mass);
			extensive[i].energy = 0.5*extensive[i].momentum*extensive[i].momentum / extensive[i].mass +
				extensive[i].mass*eos.dp2e(cells[i].density, cells[i].pressure);
			cells[i].entropy = eos.dp2s(cells[i].density, cells[i].pressure);
		}
	}
}
void hdsim::TimeAdvance2()
{
	double dt = GetTimeStep(cells_,edges_,eos_,cfl_);

	interpolation_.GetInterpolatedValues(cells_, edges_, interp_values_);
	GetRSvalues(interp_values_, rs_, rs_values_);

	vector<Extensive> old_extensive(extensives_);
	vector<double> old_edges(edges_);

	UpdateExtensives(extensives_, rs_values_, 0.5*dt);
	source_.CalcForce(edges_, cells_, time_, extensives_,0.5*dt);
	UpdateEdges(edges_, rs_values_, 0.5*dt);
	UpdateCells(extensives_, edges_, eos_, cells_,rs_values_);
	time_ += 0.5*dt;

	interpolation_.GetInterpolatedValues(cells_, edges_, interp_values_);
	GetRSvalues(interp_values_, rs_, rs_values_);

	extensives_ = old_extensive;
	edges_ = old_edges;
	UpdateExtensives(extensives_, rs_values_, dt);
	source_.CalcForce(edges_, cells_, time_, extensives_,dt);
	UpdateEdges(edges_, rs_values_, dt);
	UpdateCells(extensives_, edges_, eos_, cells_,rs_values_);
	time_ += 0.5*dt;
	++cycle_;
}

double hdsim::GetTime() const
{
	return time_;
}

vector<Primitive> const & hdsim::GetCells() const
{
	return cells_;
}

vector<double> const & hdsim::GetEdges() const
{
	return edges_;
}

size_t hdsim::GetCycle() const
{
	return cycle_;
}

void hdsim::SetTime(double t)
{
	time_ = t;
}
