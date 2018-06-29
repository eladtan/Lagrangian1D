#include "hdsim.hpp"
#include <algorithm>
#include <iostream>

hdsim::hdsim(double cfl, vector<Primitive> const& cells, vector<double> const& edges, Interpolation const& interp,
	IdealGas const& eos, ExactRS const& rs,SourceTerm const& source, Geometry const& geo, 
	BoundarySolution const* BS):cfl_(cfl),
	cells_(cells),edges_(edges),interpolation_(interp),eos_(eos),rs_(rs),time_(0),cycle_(0),
	extensives_(vector<Extensive>()),source_(source),geo_(geo), BoundarySolution_(BS)
{
	size_t N = cells.size();
	extensives_.resize(N);
	for (size_t i = 0; i < N; ++i)
	{
		double vol = geo_.GetVolume(edges,i);;
		extensives_[i].mass = cells[i].density*vol;
		extensives_[i].momentum = extensives_[i].mass*cells[i].velocity;
		extensives_[i].energy = 0.5*extensives_[i].momentum*extensives_[i].momentum / extensives_[i].mass +
			eos_.dp2e(cells[i].density, cells[i].pressure)*extensives_[i].mass;
		extensives_[i].entropy = eos_.dp2s(cells[i].density, cells[i].pressure)*extensives_[i].mass;
		cells_[i].energy = eos_.dp2e(cells_[i].density, cells_[i].pressure);
		extensives_[i].et = cells_[i].energy*extensives_[i].mass;
	}
}


hdsim::~hdsim()
{}

namespace
{
	double GetTimeStep(vector<Primitive> const& cells, vector<double> const& edges, IdealGas const& eos, double cfl,
		SourceTerm const&source,std::vector<RSsolution> const& rs_values)
	{
		double force_inverse_dt = source.GetInverseTimeStep(edges);
		double dt_1 = eos.dp2c(cells[0].density, cells[0].pressure) / (edges[1] - edges[0]);
		size_t N = cells.size();
		for (size_t i = 1; i < N-1; ++i)
		{
			double dv = -(rs_values[i + 1].velocity - rs_values[i].velocity);
			dt_1 = std::max(dt_1, std::max(dv, eos.dp2c(cells[i].density, cells[i].pressure)) / (edges[i + 1] - edges[i]));
		}
		dt_1 = std::max(dt_1, eos.dp2c(cells[N - 1].density, cells[N - 1].pressure) / (edges[N] - edges[N-1]));
		dt_1 = max(dt_1, force_inverse_dt);
		return cfl/dt_1;
	}

	void GetRSvalues(vector<pair<Primitive,Primitive> > const& interp_values, ExactRS const&rs,
		vector<RSsolution> &res)
	{
		size_t N = interp_values.size();
		res.resize(N);
		for (size_t i = 0; i < N; ++i)
			res[i] = rs.Solve(interp_values[i].first, interp_values[i].second);
	}

	void UpdateExtensives(vector<Extensive> &cells, vector<RSsolution> const& rs_values_,double dt,Geometry const& geo,
		std::vector<double> const& edges)
	{
		size_t N = cells.size();
		for (size_t i = 0; i < N; ++i)
		{
			double v = cells[i].momentum / cells[i].mass;
			double dP = (rs_values_[i + 1].pressure * geo.GetArea(edges[i+1]) - rs_values_[i].pressure 
				* geo.GetArea(edges[i]));
			cells[i].momentum -= dP*dt;
			double dE = (rs_values_[i + 1].pressure*rs_values_[i + 1].velocity* geo.GetArea(edges[i + 1])
				- rs_values_[i].pressure*rs_values_[i].velocity* geo.GetArea(edges[i]))*dt;
			cells[i].energy -= dE;
			cells[i].et -= dE - v*dP*dt;
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
		double dv = rsvalues[index + 1].velocity - rsvalues[index].velocity;
		if (dv > 0)
			return true;
		double ek = cell.velocity*cell.velocity;
		if (et < 0)
			return true;
		if (et > 0.0001*ek)
			return false;
/*		double de = rsvalues[index + 1].velocity*rsvalues[index+1].pressure - rsvalues[index].velocity*
			rsvalues[index].pressure;
		if (de<0)
			return true;*/
		if (dv*dv > 0.0001*et)
			return false;
		else
			return true;
	}

	void UpdateCells(vector<Extensive> &extensive, vector<double> const& edges, IdealGas const& eos,
		vector<Primitive> &cells,vector<RSsolution> const& rsvalues, Geometry const& geo)
	{
		size_t N = cells.size();
		for (size_t i = 0; i < N; ++i)
		{
			double vol = geo.GetVolume(edges, i);
			cells[i].density = extensive[i].mass / vol;
			cells[i].velocity = extensive[i].momentum / extensive[i].mass;
			cells[i].energy = extensive[i].et / extensive[i].mass;
			double et = cells[i].energy;
			if (ShouldUseEntropy(cells[i], rsvalues, i,et))
				cells[i].pressure = eos.sd2p(extensive[i].entropy/extensive[i].mass, cells[i].density);
			else
				cells[i].pressure = eos.de2p(cells[i].density, et);
			et = extensive[i].mass*eos.dp2e(cells[i].density, cells[i].pressure);
			extensive[i].energy = 0.5*extensive[i].momentum*extensive[i].momentum / extensive[i].mass +	et;
			extensive[i].et = et;
			cells[i].entropy = eos.dp2s(cells[i].density, cells[i].pressure);
			extensive[i].entropy = cells[i].entropy*extensive[i].mass;
		}
	}
}

void hdsim::ReCalcCells(vector<Extensive> const& extensive)
{
	size_t N = cells_.size();
	for (size_t i = 0; i < N; ++i)
	{
		double vol = geo_.GetVolume(edges_, i);
		cells_[i].density = extensive[i].mass / vol;
		cells_[i].velocity = extensive[i].momentum / extensive[i].mass;
		const double et = extensive[i].et / extensive[i].mass;
		cells_[i].pressure = eos_.de2p(cells_[i].density, et);
		cells_[i].entropy = eos_.dp2s(cells_[i].density, cells_[i].pressure);
	}
}


void hdsim::TimeAdvance()
{
	interpolation_.GetInterpolatedValues(cells_, edges_, interp_values_,time_);
	GetRSvalues(interp_values_, rs_, rs_values_);
	double dt = GetTimeStep(cells_, edges_, eos_, cfl_, source_,rs_values_);


	if (BoundarySolution_ != 0)
	{
		pair<RSsolution,RSsolution> bvalues = BoundarySolution_->GetBoundaryValues(cells_);
		if (BoundarySolution_->ShouldCalc().first)
			rs_values_[0] = bvalues.first;
		if (BoundarySolution_->ShouldCalc().second)
			rs_values_.back() = bvalues.second;
	}

	UpdateExtensives(extensives_, rs_values_, dt,geo_,edges_);
	source_.CalcForce(edges_, cells_, time_, extensives_,dt);
	UpdateEdges(edges_, rs_values_, dt);
	UpdateCells(extensives_, edges_, eos_, cells_,rs_values_,geo_);
	time_ += dt;
	++cycle_;
}

void hdsim::TimeAdvance2()
{
	interpolation_.GetInterpolatedValues(cells_, edges_, interp_values_, time_);
	GetRSvalues(interp_values_, rs_, rs_values_);

	if (BoundarySolution_ != 0)
	{
		pair<RSsolution, RSsolution> bvalues = BoundarySolution_->GetBoundaryValues(cells_);
		if (BoundarySolution_->ShouldCalc().first)
			rs_values_[0] = bvalues.first;
		if (BoundarySolution_->ShouldCalc().second)
			rs_values_.back() = bvalues.second;
	}
	double dt = GetTimeStep(cells_, edges_, eos_, cfl_, source_,rs_values_);
	if (cycle_ == 0)
		dt *= 0.005;

	vector<Extensive> old_extensive(extensives_);
	vector<double> old_edges(edges_);

	UpdateExtensives(extensives_, rs_values_, 0.5*dt,geo_,edges_);
	source_.CalcForce(edges_, cells_, time_, extensives_, 0.5*dt);
	UpdateEdges(edges_, rs_values_, 0.5*dt);
	UpdateCells(extensives_, edges_, eos_, cells_, rs_values_,geo_);
	time_ += 0.5*dt;

	interpolation_.GetInterpolatedValues(cells_, edges_, interp_values_, time_);
	GetRSvalues(interp_values_, rs_, rs_values_);
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
	UpdateExtensives(extensives_, rs_values_, dt,geo_,edges_);
	source_.CalcForce(edges_, cells_, time_, extensives_, dt);
	UpdateEdges(edges_, rs_values_, dt);
	UpdateCells(extensives_, edges_, eos_, cells_, rs_values_,geo_);
	time_ += 0.5*dt;
	++cycle_;
}


void hdsim::SetCycle(size_t cyc)
{
	cycle_ = cyc;
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
