#include "hdsim.hpp"
#include <algorithm>
#include <iostream>

hdsim::hdsim(double cfl, vector<Primitive> const& cells, vector<double> const& edges, MinMod const& interp,
	IdealGas const& eos, ExactRS const& rs,SourceTerm const& source,BoundarySolution const* BS):cfl_(cfl),
	cells_(cells),edges_(edges),interpolation_(interp),eos_(eos),rs_(rs),time_(0),cycle_(0),
	extensives_(vector<Extensive>()),source_(source),BoundarySolution_(BS)
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
		extensives_[i].entropy = eos_.dp2s(cells[i].density, cells[i].pressure)*extensives_[i].mass;
	}
}


hdsim::~hdsim()
{}

namespace
{
	double GetTimeStep(vector<Primitive> const& cells, vector<double> const& edges, IdealGas const& eos, double cfl,
		SourceTerm const&source)
	{
		double force_inverse_dt = source.GetInverseTimeStep(edges);
		double dt = (edges[1] - edges[0]) / eos.dp2c(cells[0].density, cells[0].pressure);
		size_t N = cells.size();
		for (size_t i = 1; i < N;++i)
			dt = min(dt,(edges[i+1] - edges[i]) / eos.dp2c(cells[i].density, cells[i].pressure));
		dt = max(1.0 / dt, force_inverse_dt);
		return cfl/dt;
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

	bool ShouldUseEntropy(Primitive const& cell, vector<RSsolution> const& rsvalues, size_t index,double et)
	{
		double ek = cell.velocity*cell.velocity;
		if (et < 0)
			return true;
		if (et > 0.01*ek)
			return false;
		double dv = rsvalues[index + 1].velocity - rsvalues[index].velocity;
		if (dv > 0)
			return true;
		double de = rsvalues[index + 1].velocity*rsvalues[index+1].pressure - rsvalues[index].velocity*
			rsvalues[index].pressure;
		if (de<0)
			return true;
		if (dv*dv > 0.3*et)
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
			const double et = (extensive[i].energy - 0.5*extensive[i].momentum*extensive[i].momentum
				/ extensive[i].mass) / extensive[i].mass;
			if (ShouldUseEntropy(cells[i], rsvalues, i,et))
				cells[i].pressure = eos.sd2p(extensive[i].entropy/extensive[i].mass, cells[i].density);
			else
				cells[i].pressure = eos.de2p(cells[i].density, et);
			extensive[i].energy = 0.5*extensive[i].momentum*extensive[i].momentum / extensive[i].mass +
				extensive[i].mass*eos.dp2e(cells[i].density, cells[i].pressure);
			cells[i].entropy = eos.dp2s(cells[i].density, cells[i].pressure);
			extensive[i].entropy = cells[i].entropy*extensive[i].mass;
		}
	}
}
void hdsim::TimeAdvance2()
{
	double dt = GetTimeStep(cells_,edges_,eos_,cfl_,source_);

	interpolation_.GetInterpolatedValues(cells_, edges_, interp_values_,time_);
	GetRSvalues(interp_values_, rs_, rs_values_);

	if (BoundarySolution_ != 0)
	{
		pair<RSsolution,RSsolution> bvalues = BoundarySolution_->GetBoundaryValues(cells_);
		if (BoundarySolution_->ShouldCalc().first)
			rs_values_[0] = bvalues.first;
		if (BoundarySolution_->ShouldCalc().second)
			rs_values_.back() = bvalues.second;
	}

	vector<Extensive> old_extensive(extensives_);
	vector<double> old_edges(edges_);

	UpdateExtensives(extensives_, rs_values_, 0.5*dt);
	source_.CalcForce(edges_, cells_, time_, extensives_,0.5*dt);
	UpdateEdges(edges_, rs_values_, 0.5*dt);
	UpdateCells(extensives_, edges_, eos_, cells_,rs_values_);
	time_ += 0.5*dt;

	interpolation_.GetInterpolatedValues(cells_, edges_, interp_values_,time_);
	GetRSvalues(interp_values_, rs_, rs_values_);
//	std::cout << BoundarySolution_->ShouldCalc().second << std::endl;
	if (BoundarySolution_ != 0)
	{
		pair<RSsolution, RSsolution> bvalues = BoundarySolution_->GetBoundaryValues(cells_);
		if (BoundarySolution_->ShouldCalc().first)
			rs_values_[0] = bvalues.first;
		if (BoundarySolution_->ShouldCalc().second)
		{
			rs_values_.back() = bvalues.second;
			//std::cout << rs_values_.back().velocity << std::endl;
		}
	}

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

vector<Primitive>& hdsim::GetCells()
{
	return cells_;
}

vector<Extensive> const & hdsim::GetExtensives() const
{
	return extensives_;
}

vector<Extensive> & hdsim::GetExtensives() 
{
	return extensives_;
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
