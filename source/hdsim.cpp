#include "hdsim.hpp"
#include <algorithm>
#include <iostream>
#ifdef RICH_MPI
#include "mpi_comm.hpp"
#endif
namespace
{
	template <class T> vector<T> unique(vector<T> const& v)
	{
		std::size_t n = v.size();
		vector<T> res;
		res.reserve(n);
		if (n == 0)
			return res;
		res.push_back(v[0]);
		for (typename vector<T>::const_iterator it = v.begin() + 1; it != v.end(); ++it)
			if (*it == *(it - 1))
				continue;
			else
				res.push_back(*it);
		return res;
	}


	template <class T> void RemoveVector
	(vector<T> &v, vector<size_t> &indeces)
	{
		if (indeces.empty())
			return;
		sort(indeces.begin(), indeces.end());
		vector<T> result;
		result.reserve(v.size() - indeces.size());
		int counter = 0;
		for (std::size_t i = 0; i < static_cast<std::size_t>(indeces.back()); ++i)
		{
			if (std::size_t(indeces[std::size_t(counter)]) == i)
				++counter;
			else
				result.push_back(v[i]);
		}
		for (std::size_t i = static_cast<std::size_t>(indeces.back()) + 1; i < v.size(); ++i)
			result.push_back(v[i]);
		v = result;
	}
}

hdsim::hdsim(double cfl, vector<Primitive> const& cells, vector<double> const& edges, Interpolation const& interp,
	IdealGas const& eos, RiemannSolver const& rs,SourceTerm const& source, Geometry const& geo,
	const double AMR_ratio, BoundarySolution const* BS):cfl_(cfl),
	cells_(cells),edges_(edges),interpolation_(interp),eos_(eos),rs_(rs),time_(0),cycle_(0),TotalEcool_(0),
	extensives_(vector<Extensive>()),source_(source),geo_(geo), AMR_ratio_(AMR_ratio), BoundarySolution_(BS),dt_suggest_(-1)
{
#ifdef RICH_MPI
	size_t Ntotal = cells.size();
	int rank = 0,ws=0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	size_t index_lower = rank * Ntotal / ws;
	size_t index_upper = (rank + 1)*Ntotal / ws;
	std::vector<Primitive> cells_temp(cells.begin() + index_lower, 
		cells.begin() + index_upper);
	std::vector<Primitive> cells_temp2(cells.begin(), cells.begin() + 1);
	cells_ = cells_temp;
	std::vector<double> edge_temp(edges.begin() + index_lower, edges.begin() +
		index_upper + 1);
	edges_ = edge_temp;
#endif
	size_t N = cells_.size();
	extensives_.resize(N);
	for (size_t i = 0; i < N; ++i)
	{
		double vol = geo_.GetVolume(edges_,i);;
		extensives_[i].mass = cells_[i].density*vol;
		extensives_[i].momentum = extensives_[i].mass*cells_[i].velocity;
		extensives_[i].energy = 0.5*extensives_[i].momentum*extensives_[i].momentum / extensives_[i].mass +
			eos_.dp2e(cells_[i].density, cells_[i].pressure)*extensives_[i].mass;
		extensives_[i].entropy = eos_.dp2s(cells_[i].density, cells_[i].pressure)*extensives_[i].mass;
		cells_[i].energy = eos_.dp2e(cells_[i].density, cells_[i].pressure);
		cells_[i].LastCool = time_;
		extensives_[i].et = cells_[i].energy*extensives_[i].mass;
	}
}


hdsim::~hdsim()
{}

namespace
{
	double GetTimeStep(vector<Primitive> const& cells, vector<double> const& edges, IdealGas const& eos, double cfl,
		SourceTerm const&source,std::vector<RSsolution> const& rs_values,double &dt_suggest)
	{
		double force_inverse_dt = source.GetInverseTimeStep(edges);
		double dt_1 = std::max(rs_values[0].velocity - rs_values[1].velocity,eos.dp2c(cells[0].density, cells[0].pressure)) / (edges[1] - edges[0]);
		size_t N = cells.size();
		for (size_t i = 1; i < N; ++i)
		{
			double dv = -(rs_values[i + 1].velocity - rs_values[i].velocity);
			dt_1 = std::max(dt_1, std::max(dv, eos.dp2c(cells[i].density, cells[i].pressure)) / (edges[i + 1] - edges[i]));
		}
		//dt_1 = std::max(dt_1, eos.dp2c(cells[N - 1].density, cells[N - 1].pressure) / (edges[N] - edges[N-1]));
		dt_1 = max(dt_1, force_inverse_dt);
		if (dt_suggest > 0)
		{
			dt_1 = std::max(dt_1, 1.0 / dt_suggest);
			dt_suggest = -1;
		}
#ifdef RICH_MPI
		MPI_Allreduce(MPI_IN_PLACE, &dt_1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
		return cfl/dt_1;
	}

	void GetRSvalues(vector<pair<Primitive,Primitive> > const& interp_values, RiemannSolver const&rs,
		vector<RSsolution> &res)
	{
		size_t N = interp_values.size();
		res.resize(N);
		for (int i = 0; i < N; ++i)
			res[i] = rs.Solve(interp_values[i].first, interp_values[i].second);
	}

	void UpdateExtensives(vector<Extensive> &cells, vector<RSsolution> const& rs_values_,double dt,Geometry const& geo,
		std::vector<double> const& edges)
	{
		size_t N = cells.size();
		for (size_t i = 0; i < N; ++i)
		{
			double v = cells[i].momentum / cells[i].mass;
			double dP = -(rs_values_[i + 1].pressure * geo.GetArea(edges[i + 1]) - 
				rs_values_[i].pressure * geo.GetArea(edges[i]));
			double oldEk = cells[i].momentum*cells[i].momentum / cells[i].mass;
			cells[i].momentum += dP * dt;
			double dE = -(rs_values_[i + 1].pressure*rs_values_[i+1].velocity* geo.GetArea(edges[i + 1]) - 
				rs_values_[i].pressure*rs_values_[i].velocity* geo.GetArea(edges[i]))*dt;
			cells[i].energy += dE;
			double newEk = cells[i].momentum*cells[i].momentum / cells[i].mass;
			cells[i].et += dE - 0.5*(newEk - oldEk);
		}
	}

	void UpdateEdges(vector<double> &edges, vector<RSsolution> const& rs_values_,double dt)
	{
		size_t N = edges.size();
		for (size_t i = 0; i < N; ++i)
			edges[i] += rs_values_[i].velocity*dt;
	}

	bool ShouldUseEntropy(Primitive const& cell, vector<RSsolution> const& rsvalues, size_t index,double et,double EntropyEt)
	{
		double dv = rsvalues[index + 1].velocity - rsvalues[index].velocity;
	//	if (dv > 0)
	//		return true;
		//if (et < EntropyEt)
		//	return true;
		double ek = cell.velocity*cell.velocity;
		if (et < 0)
			return true;
		if (et > 0.001*ek)
			return false;
		if (dv*dv > 0.001*et && dv <0)
			return false;
		else
			return true;
	}

	void UpdateCells(vector<Extensive> &extensive, vector<double> const& edges, IdealGas const& eos,
		vector<Primitive> &cells,vector<RSsolution> const& rsvalues, Geometry const& geo)
	{
		size_t N = cells.size();
		for (int i = 0; i < N; ++i)
		{
			double vol = geo.GetVolume(edges, i);
			cells[i].density = extensive[i].mass / vol;
			cells[i].velocity = extensive[i].momentum / extensive[i].mass;
			cells[i].energy = extensive[i].et / extensive[i].mass;
			double et = cells[i].energy;
			double et2 = extensive[i].energy / extensive[i].mass - 0.5 * cells[i].velocity*cells[i].velocity;
			if (ShouldUseEntropy(cells[i], rsvalues, i,et,eos.dp2e(cells[i].density,eos.sd2p(cells[i].entropy,cells[i].density))))
				cells[i].pressure = eos.sd2p(extensive[i].entropy/extensive[i].mass, cells[i].density);
			else
				cells[i].pressure = eos.de2p(cells[i].density, et);
			cells[i].pressure = std::max(cells[i].pressure, cells[i].density*cells[i].velocity*cells[i].velocity*0.0001);
			et = extensive[i].mass*eos.dp2e(cells[i].density, cells[i].pressure);
			extensive[i].energy = 0.5*extensive[i].momentum*extensive[i].momentum / extensive[i].mass +	et;
			extensive[i].et = et;
			cells[i].entropy = eos.dp2s(cells[i].density, cells[i].pressure);
			extensive[i].entropy = cells[i].entropy*extensive[i].mass;
		}
	}
}

void hdsim::ReCalcCells(vector<Extensive> &extensive)
{
	size_t N = cells_.size();
	for (size_t i = 0; i < N; ++i)
	{
		double vol = geo_.GetVolume(edges_, i);
		cells_[i].density = extensive[i].mass / vol;
		cells_[i].velocity = extensive[i].momentum / extensive[i].mass;
		const double et = extensive[i].et / extensive[i].mass;
		cells_[i].pressure = eos_.de2p(cells_[i].density, et);
		cells_[i].energy = et;
		cells_[i].entropy = eos_.dp2s(cells_[i].density, cells_[i].pressure);
		extensive[i].entropy = cells_[i].entropy*extensive[i].mass;
	}
}

void hdsim::ReCalcExtensives(vector<Primitive> const& cells)
{
	size_t N = cells.size();
	for (size_t i = 0; i < N; ++i)
	{
		double vol = geo_.GetVolume(edges_, i);
		extensives_[i].mass = vol*cells[i].density;
		extensives_[i].momentum = extensives_[i].mass*cells[i].velocity;
		extensives_[i].et = extensives_[i].mass*cells[i].energy;
		extensives_[i].entropy = extensives_[i].mass*cells[i].entropy;
		extensives_[i].energy = extensives_[i].et + 0.5*extensives_[i].momentum*extensives_[i].momentum / extensives_[i].mass;
	}
}

void hdsim::SuggestTimeStep(double dt)
{
	dt_suggest_ = dt;
}

void hdsim::TimeAdvance()
{
	interpolation_.GetInterpolatedValues(cells_, edges_, interp_values_,time_);
	GetRSvalues(interp_values_, rs_, rs_values_);
	double dt = GetTimeStep(cells_, edges_, eos_, cfl_, source_,rs_values_,dt_suggest_);

	if (BoundarySolution_ != 0)
	{
#ifdef RICH_MPI
		int rank = 0, ws = 0;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &ws);
		if (rank == 0 || rank == (ws - 1))
		{
#endif
			pair<RSsolution, RSsolution> bvalues = BoundarySolution_->GetBoundaryValues(cells_);
#ifdef RICH_MPI
			if (rank == 0)
#endif
				if (BoundarySolution_->ShouldCalc().first)
					rs_values_[0] = bvalues.first;
#ifdef RICH_MPI
			if (rank == ws - 1)
#endif
				if (BoundarySolution_->ShouldCalc().second)
					rs_values_.back() = bvalues.second;
#ifdef RICH_MPI
		}
#endif
	}
	UpdateExtensives(extensives_, rs_values_, dt,geo_,edges_);
	source_.CalcForce(edges_, cells_, time_, extensives_,dt);
	UpdateEdges(edges_, rs_values_, dt);
	UpdateCells(extensives_, edges_, eos_, cells_,rs_values_,geo_);
	time_ += dt;
	++cycle_;
	AMR();
}

void hdsim::TimeAdvance2()
{
	interpolation_.GetInterpolatedValues(cells_, edges_, interp_values_, time_);
	GetRSvalues(interp_values_, rs_, rs_values_);

	if (BoundarySolution_ != 0)
	{
#ifdef RICH_MPI
		int rank = 0,ws=0;
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		MPI_Comm_size(MPI_COMM_WORLD,&ws);
		if (rank == 0 || rank == (ws - 1))
		{
#endif
			pair<RSsolution, RSsolution> bvalues = BoundarySolution_->GetBoundaryValues(cells_);
#ifdef RICH_MPI
			if (rank == 0)
#endif
				if (BoundarySolution_->ShouldCalc().first)
					rs_values_[0] = bvalues.first;
#ifdef RICH_MPI
			if (rank == ws - 1)
#endif
				if (BoundarySolution_->ShouldCalc().second)
					rs_values_.back() = bvalues.second;
#ifdef RICH_MPI
		}
#endif
	}
	double dt = GetTimeStep(cells_, edges_, eos_, cfl_, source_,rs_values_,dt_suggest_);
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
#ifdef RICH_MPI
		int rank = 0, ws = 0;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &ws);
		if (rank == 0 || rank == (ws - 1))
		{
#endif
			pair<RSsolution, RSsolution> bvalues = BoundarySolution_->GetBoundaryValues(cells_);
#ifdef RICH_MPI
			if (rank == 0)
#endif
				if (BoundarySolution_->ShouldCalc().first)
					rs_values_[0] = bvalues.first;
#ifdef RICH_MPI
			if (rank == ws - 1)
#endif
				if (BoundarySolution_->ShouldCalc().second)
					rs_values_.back() = bvalues.second;
#ifdef RICH_MPI
		}
#endif
	}
	extensives_ = old_extensive;
	edges_ = old_edges;
	UpdateExtensives(extensives_, rs_values_, dt,geo_,edges_);
	source_.CalcForce(edges_, cells_, time_, extensives_, dt);
	UpdateEdges(edges_, rs_values_, dt);
#ifdef RICH_MPI
	RedistributeExtensives(extensives_,edges_,cells_,rs_values_);
#endif
	UpdateCells(extensives_, edges_, eos_, cells_, rs_values_,geo_);
	time_ += 0.5*dt;
	++cycle_;
	AMR();
}

void hdsim::AMR(void)
{
	double pratio = 1.04;
	double dratio = 1.1;
	size_t N = cells_.size();
	std::vector<size_t> edge_remove;
	for (size_t i = 3; i < N - 3; ++i)
	{
		// was left cell removed?
		if (edge_remove.size() > 0 && (edge_remove.back() == i || edge_remove.back() == i - 1))
			continue;
		double dx = edges_[i + 1] - edges_[i];
		// Are we too small?
		if(dx< AMR_ratio_*edges_[i])
		{
			// Are we smooth?
			bool smooth_left_left_left = (cells_[i - 2].density < cells_[i - 3].density * dratio) && (cells_[i - 2].density * dratio > cells_[i - 3].density)
				&& (cells_[i - 2].pressure < cells_[i - 3].pressure*pratio) && (cells_[i - 2].pressure*pratio > cells_[i - 3].pressure);
			bool smooth_left_left = (cells_[i - 1].density < cells_[i - 2].density * dratio) && (cells_[i - 1].density * dratio > cells_[i - 2].density)
				&& (cells_[i - 1].pressure < cells_[i - 2].pressure*pratio) && (cells_[i - 1].pressure*pratio > cells_[i - 2].pressure);
			bool smooth_left = (cells_[i].density < cells_[i - 1].density * dratio) && (cells_[i].density * dratio > cells_[i - 1].density)
				&& (cells_[i].pressure < cells_[i - 1].pressure*pratio) && (cells_[i].pressure*pratio > cells_[i - 1].pressure);
			bool smooth_right = (cells_[i].density < cells_[i + 1].density * dratio) && (cells_[i].density * dratio > cells_[i + 1].density)
				&& (cells_[i].pressure < cells_[i + 1].pressure*pratio) && (cells_[i].pressure*pratio > cells_[i + 1].pressure);
			bool smooth_right_right = (cells_[i + 1].density < cells_[i + 2].density * dratio) && (cells_[i + 1].density * dratio > cells_[i + 2].density)
				&& (cells_[i + 1].pressure < cells_[i + 2].pressure*pratio) && (cells_[i + 1].pressure*pratio > cells_[i + 2].pressure);
			bool smooth_right_right_right = (cells_[i + 2].density < cells_[i + 3].density * dratio) && (cells_[i + 2].density * dratio > cells_[i + 3].density)
				&& (cells_[i + 2].pressure < cells_[i + 3].pressure*pratio) && (cells_[i + 2].pressure*pratio > cells_[i + 3].pressure);
			if (smooth_left && smooth_right && smooth_left_left && smooth_right_right && smooth_left_left_left && smooth_right_right_right)
			{
				if((edges_[i+2]-edges_[i+1])>(edges_[i] - edges_[i - 1]))
					edge_remove.push_back(i);
				else
					edge_remove.push_back(i + 1);
			}
		}
	}
	edge_remove = unique(edge_remove);
	size_t Nremove = edge_remove.size();
	for (size_t i = 0; i < Nremove; ++i)
		extensives_[edge_remove[i] - 1] += extensives_[edge_remove[i]];
	// Remove old cells
	if (Nremove > 0)
	{
		RemoveVector(extensives_, edge_remove);
		RemoveVector(cells_, edge_remove);
		RemoveVector(edges_, edge_remove);
		ReCalcCells(extensives_);
	}
}


void hdsim::SetCycle(size_t cyc)
{
	cycle_ = cyc;
}

double hdsim::GetTime() const
{
	return time_;
}

double hdsim::GetEcool() const
{
	return TotalEcool_;
}

double& hdsim::GetEcool()
{
	return TotalEcool_;
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


vector<double>& hdsim::GetEdges() 
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

void hdsim::SetEcool(double E)
{
	TotalEcool_ = E;
}
