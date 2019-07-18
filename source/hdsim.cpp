#include "hdsim.hpp"
#include <algorithm>
#include <iostream>
#ifdef RICH_MPI
#include "mpi_comm.hpp"
#endif

namespace
{
	Extensive CalcExtensive(Geometry const& geo,std::vector<double> const&
		edges,Primitive const& cell,size_t index)
	{
		Extensive res;
		double vol = geo.GetVolume(edges, index);
		res.mass = vol * cell.density;
		res.momentum = res.mass*cell.velocity;
		res.et = res.mass*cell.energy;
		res.entropy = res.mass*cell.entropy;
		res.energy = res.et + 0.5*res.momentum*res.momentum / res.mass;
		return res;
	}

	Primitive CalcPrimitive(double vol, Extensive const& extensive,
		IdealGas const& eos)
	{
		Primitive cell;
		cell.density = extensive.mass / vol;
		cell.velocity = extensive.momentum / extensive.mass;
		cell.energy = extensive.et / extensive.mass;
		cell.pressure = eos.de2p(cell.density, cell.energy);
		cell.entropy = eos.dp2s(cell.density, cell.pressure);
		return cell;
	}

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
	IdealGas const& eos, RiemannSolver const& rs, SourceTerm const& source, Geometry const& geo,
	const double AMR_ratio, BoundarySolution const* BS) :cfl_(cfl),
	cells_(cells), edges_(edges), interpolation_(interp), eos_(eos), rs_(rs), time_(0), cycle_(0), TotalEcool_(0),
	extensives_(vector<Extensive>()), source_(source), geo_(geo), AMR_ratio_(AMR_ratio), BoundarySolution_(BS), dt_suggest_(-1)
{
#ifdef RICH_MPI
	size_t Ntotal = cells.size();
	int rank = 0, ws = 0;
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
		double vol = geo_.GetVolume(edges_, i);;
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

	double GetTimeStep(vector<Primitive> const& cells, vector<double> const& edges, IdealGas const& eos, double cfl,
		SourceTerm const&source, std::vector<RSsolution> const& rs_values, double &dt_suggest)
	{
		double force_inverse_dt = source.GetInverseTimeStep(edges);
		double dt_1 = std::max(3 * (rs_values[0].velocity - rs_values[1].velocity), eos.dp2c(cells[0].density, cells[0].pressure)) / (edges[1] - edges[0]);
		size_t N = cells.size();
		for (size_t i = 1; i < N; ++i)
		{
			double dv = -(rs_values[i + 1].velocity - rs_values[i].velocity);
			dt_1 = std::max(dt_1, std::max(3 * dv, eos.dp2c(cells[i].density, cells[i].pressure)) / (edges[i + 1] - edges[i]));
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
		return cfl / dt_1;
	}

	void GetRSvalues(vector<pair<Primitive, Primitive> > const& interp_values, RiemannSolver const&rs,
		vector<RSsolution> &res)
	{
		size_t N = interp_values.size();
		res.resize(N);
		for (int i = 0; i < N; ++i)
			res[i] = rs.Solve(interp_values[i].first, interp_values[i].second);
	}

	void UpdateExtensives(vector<Extensive> &cells, vector<RSsolution> const& rs_values_, double dt, Geometry const& geo,
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
			double dE = -(rs_values_[i + 1].pressure*rs_values_[i + 1].velocity* geo.GetArea(edges[i + 1]) -
				rs_values_[i].pressure*rs_values_[i].velocity* geo.GetArea(edges[i]))*dt;
			cells[i].energy += dE;
			double newEk = cells[i].momentum*cells[i].momentum / cells[i].mass;
			cells[i].et += dE - 0.5*(newEk - oldEk);
		}
	}

	void UpdateEdges(vector<double> &edges, vector<RSsolution> const& rs_values_, double dt)
	{
		size_t N = edges.size();
		for (size_t i = 0; i < N; ++i)
			edges[i] += rs_values_[i].velocity*dt;
	}

	bool ShouldUseEntropy(Primitive const& cell, vector<RSsolution> const& rsvalues, size_t index, double et, double EntropyEt)
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
		if (dv*dv > 0.001*et && dv < 0)
			return false;
		else
			return true;
	}

	void UpdateCells(vector<Extensive> &extensive, vector<double> const& edges, IdealGas const& eos,
		vector<Primitive> &cells, vector<RSsolution> const& rsvalues, Geometry const& geo)
	{
		size_t N = cells.size();
		for (int i = 0; i < N; ++i)
		{
			double vol = geo.GetVolume(edges, i);
			if (vol < 0)
			{
				std::cout << "Negative vol = " << vol << " x = " << edges[i] << " index = " << i << " Rs p/v left " << rsvalues[i].pressure << "," << rsvalues[i].velocity <<
					" Rs p/v right " << rsvalues[i + 1].pressure << "," << rsvalues[i + 1].velocity << std::endl;
				std::cout << " left rho " << cells[i].density << " pressure " << cells[i].pressure << " v " << cells[i].velocity
					<< " right rho " << cells[i + 1].density << " pressure " << cells[i + 1].pressure << " v " << cells[i + 1].velocity << std::endl;
				throw;
			}
			cells[i].density = extensive[i].mass / vol;
			cells[i].velocity = extensive[i].momentum / extensive[i].mass;
			cells[i].energy = extensive[i].et / extensive[i].mass;
			double et = cells[i].energy;
			double et2 = extensive[i].energy / extensive[i].mass - 0.5 * cells[i].velocity*cells[i].velocity;
			if (ShouldUseEntropy(cells[i], rsvalues, i, et, eos.dp2e(cells[i].density, eos.sd2p(cells[i].entropy, cells[i].density))))
				cells[i].pressure = eos.sd2p(extensive[i].entropy / extensive[i].mass, cells[i].density);
			else
				cells[i].pressure = eos.de2p(cells[i].density, et);
			//			cells[i].pressure = std::max(cells[i].pressure, cells[i].density*cells[i].velocity*cells[i].velocity*0.0001);
			et = extensive[i].mass*eos.dp2e(cells[i].density, cells[i].pressure);
			extensive[i].energy = 0.5*extensive[i].momentum*extensive[i].momentum / extensive[i].mass + et;
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
		extensives_[i].mass = vol * cells[i].density;
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
#ifdef RICH_MPI
	std::array<Primitive, NGHOSTCELLS * 2> ghost_cells = SendRecvPrimitive(cells_);
	std::array<double, NGHOSTCELLS * 2> ghost_edges = SendRecvEdges(edges_);
#endif
	interpolation_.GetInterpolatedValues(cells_, edges_, interp_values_, time_
#ifdef RICH_MPI
		, ghost_cells, ghost_edges
#endif
	);
	GetRSvalues(interp_values_, rs_, rs_values_);
	double dt = GetTimeStep(cells_, edges_, eos_, cfl_, source_, rs_values_, dt_suggest_);

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
	UpdateExtensives(extensives_, rs_values_, dt, geo_, edges_);
	source_.CalcForce(edges_, cells_, time_, extensives_, dt);
	UpdateEdges(edges_, rs_values_, dt);
	UpdateCells(extensives_, edges_, eos_, cells_, rs_values_, geo_);
	time_ += dt;
	++cycle_;
	AMR(
#ifdef RICH_MPI
		ghost_cells, ghost_edges
#endif
	);
}

void hdsim::TimeAdvance2()
{
#ifdef RICH_MPI
	std::array<Primitive, NGHOSTCELLS*2> ghost_cells = SendRecvPrimitive(cells_);
	std::array<double, NGHOSTCELLS*2> ghost_edges = SendRecvEdges(edges_);
#endif
	interpolation_.GetInterpolatedValues(cells_, edges_, interp_values_, time_
#ifdef RICH_MPI
		,ghost_cells,ghost_edges
#endif
	);
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
	double dt = GetTimeStep(cells_, edges_, eos_, cfl_, source_, rs_values_, dt_suggest_);
	if (cycle_ < 10)
		dt *= 0.002;

	vector<Extensive> old_extensive(extensives_);
	vector<double> old_edges(edges_);

	UpdateExtensives(extensives_, rs_values_, 0.5*dt, geo_, edges_);
	source_.CalcForce(edges_, cells_, time_, extensives_, 0.5*dt);
	UpdateEdges(edges_, rs_values_, 0.5*dt);
	UpdateCells(extensives_, edges_, eos_, cells_, rs_values_, geo_);
	time_ += 0.5*dt;

#ifdef RICH_MPI
	ghost_cells = SendRecvPrimitive(cells_);
	ghost_edges = SendRecvEdges(edges_);
#endif

	interpolation_.GetInterpolatedValues(cells_, edges_, interp_values_, time_
#ifdef RICH_MPI
		, ghost_cells, ghost_edges
#endif
	);
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
	UpdateExtensives(extensives_, rs_values_, dt, geo_, edges_);
	source_.CalcForce(edges_, cells_, time_, extensives_, dt);
	UpdateEdges(edges_, rs_values_, dt);
#ifdef RICH_MPI
	RedistributeExtensives(extensives_, edges_, cells_, rs_values_);
#endif
	UpdateCells(extensives_, edges_, eos_, cells_, rs_values_, geo_);
	time_ += 0.5*dt;
	++cycle_;
	AMR(
#ifdef RICH_MPI
		ghost_cells,ghost_edges
#endif
	);
}

void hdsim::AMR(
#ifdef RICH_MPI
	std::array<Primitive, NGHOSTCELLS * 2> ghost_cells ,
std::array<double, NGHOSTCELLS * 2> ghost_edges
#endif
)
{
	size_t Nlevels = 6;
#ifdef RICH_MPI
	Nlevels = NGHOSTCELLS;
#endif
	double min_size = std::pow(2.0,-static_cast<double>(Nlevels));
	double pratio = 1.1;
	double dratio = 1.1;
	size_t N = cells_.size();
	// remove smooth small cells
	std::vector<size_t> edge_remove;
	for (size_t i = Nlevels; i < N - Nlevels; ++i)
	{
		// was left cell removed?
		if (edge_remove.size() > 0 && (edge_remove.back() == i || edge_remove.back() == i - 1))
			continue;
		double dx = edges_[i + 1] - edges_[i];
		// Are we too small?
		if (dx < AMR_ratio_*edges_[i] * min_size)
		{
			double pratio2 = 1.1*std::pow((AMR_ratio_*edges_[i] * min_size) / dx, 0.9);;
			double Tratio2 = 1.2*std::pow((AMR_ratio_*edges_[i] * min_size) / dx, 0.1);;
			double dratio2 = 1.2*std::pow((AMR_ratio_*edges_[i] * min_size) / dx, 1.1);
			bool smooth_left = (cells_[i].density < cells_[i - 1].density * dratio2) && (cells_[i].density * dratio2 > cells_[i - 1].density)
				&& (cells_[i].pressure < cells_[i - 1].pressure*pratio2) && (cells_[i].pressure*pratio > cells_[i - 1].pressure);
			bool smooth_right = (cells_[i].density < cells_[i + 1].density * dratio2) && (cells_[i].density * dratio2 > cells_[i + 1].density)
				&& (cells_[i].pressure < cells_[i + 1].pressure*pratio2) && (cells_[i].pressure*pratio2 > cells_[i + 1].pressure);
			bool Tratio_l = (cells_[i].density*cells_[i - 1].pressure) < (Tratio2*cells_[i - 1].density*cells_[i].pressure)
				&& (Tratio2*cells_[i].density*cells_[i - 1].pressure) > (cells_[i - 1].density*cells_[i].pressure);
			bool Tratio_r = (cells_[i].density*cells_[i + 1].pressure) < (Tratio2*cells_[i + 1].density*cells_[i].pressure)
				&& (Tratio2*cells_[i].density*cells_[i + 1].pressure) > (cells_[i + 1].density*cells_[i].pressure);
			if ((smooth_left && smooth_right && Tratio_l && Tratio_r) || (dx < AMR_ratio_*edges_[i] * min_size*0.2))
			{
				if ((edges_[i + 2] - edges_[i + 1]) > (edges_[i] - edges_[i - 1]))
					edge_remove.push_back(i);
				else
					edge_remove.push_back(i + 1);
			}
			else
			{
				if (smooth_right && Tratio_r && ((edges_[i + 2] - edges_[i + 1]) <  1.5*AMR_ratio_*edges_[i + 1]))
					edge_remove.push_back(i + 1);
				if (smooth_left && Tratio_l && ((edges_[i] - edges_[i - 1]) <  1.5*AMR_ratio_*edges_[i - 1]))
					edge_remove.push_back(i);
			}
		}
		else
		{
			if (dx < AMR_ratio_*edges_[i]*0.6)
			{
				double Tratio2 = 1.2;
				bool smooth = true;
				for (size_t j = 0; j < (2 * Nlevels - 1); ++j)
				{
					if (!((cells_[i - Nlevels + j + 1].density < cells_[i - Nlevels + j].density * dratio) && (cells_[i - Nlevels + j + 1].density * dratio > cells_[i - Nlevels + j].density)
						&& (cells_[i - Nlevels + j + 1].pressure < cells_[i - Nlevels + j].pressure*pratio) && (cells_[i - Nlevels + j + 1].pressure*pratio > cells_[i - Nlevels + j].pressure)))
					{
						smooth = false;
						break;
					}
				}
				if (smooth)
				{
					if ((edges_[i + 2] - edges_[i + 1]) > (edges_[i] - edges_[i - 1]))
						edge_remove.push_back(i);
					else
						edge_remove.push_back(i + 1);
				}
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
	// split cells ahead of non-smooth
	N = cells_.size();
	std::vector<size_t> edge_split;
	size_t Nstart = Nlevels;
	size_t Nend = N - Nlevels;
	size_t shift = 0;
	std::vector<Primitive> temp_cells;
	std::vector<double> temp_edges;
#ifdef RICH_MPI
	int rank = 0, ws = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	temp_cells.resize(N + NGHOSTCELLS * 2);
	std::copy(cells_.begin(), cells_.end(), temp_cells.begin() + NGHOSTCELLS);
	std::copy(ghost_cells.begin(), ghost_cells.begin() + NGHOSTCELLS, temp_cells.begin());
	std::copy(ghost_cells.begin() + NGHOSTCELLS, ghost_cells.end(), temp_cells.begin() + NGHOSTCELLS + N);
	temp_edges.resize(N + NGHOSTCELLS * 2+1);
	std::copy(edges_.begin(), edges_.end(), temp_edges.begin() + NGHOSTCELLS);
	std::copy(ghost_edges.begin(), ghost_edges.begin() + NGHOSTCELLS, temp_edges.begin());
	std::copy(ghost_edges.begin() + NGHOSTCELLS, ghost_edges.end(), temp_edges.begin() + NGHOSTCELLS + N+1);
	if(rank<(ws-1))
		Nend = N;
	if (rank > 0)
		Nstart = 0;
	shift = NGHOSTCELLS;
#else
	temp_cells = cells_;
	temp_edges = edges_;
#endif
	for (size_t i = Nstart; i < Nend; ++i)
	{
		double dx = temp_edges[i + 1 + shift] - temp_edges[i+shift];
		if (dx > AMR_ratio_*temp_edges[i + shift]*min_size*2)
		{
			double Tratio2 = 1.2;
			// Are we smooth?
			bool smooth = false;
			for (size_t j = 0; j < (Nlevels-1); ++j)
			{
				bool local_smooth = (temp_cells[i - Nlevels + 1 + j + shift].density < temp_cells[i - Nlevels + j + shift].density
					* dratio) && (temp_cells[i - Nlevels + 1 + j + shift].density * dratio >
						temp_cells[i - Nlevels + j + shift].density)
					&& (temp_cells[i - Nlevels + 1 + j + shift].pressure < temp_cells[i - Nlevels + j + shift].pressure
						*pratio) && (temp_cells[i - Nlevels + 1 + j + shift].pressure*pratio >
							temp_cells[i - Nlevels + j + shift].pressure);
				double l_reduce = std::pow(2.0, -static_cast<double>(j));
				if (!local_smooth && dx > AMR_ratio_*edges_[i] * 1.3*l_reduce)
				{
					smooth = true;
					break;
				}
			}
			if (!smooth)
			{
				for (size_t j = 1; j < (Nlevels - 1); ++j)
				{
					bool local_smooth = (temp_cells[i + 1 + j + shift].density < temp_cells[i + j + shift].density
						* dratio) && (temp_cells[i + 1 + j + shift].density * dratio >
							temp_cells[i + j + shift].density)
						&& (temp_cells[i + 1 + j + shift].pressure < temp_cells[i + j + shift].pressure
							*pratio) && (temp_cells[i + 1 + j + shift].pressure*pratio >
								temp_cells[i + j + shift].pressure);
					double l_reduce = std::pow(2.0, -static_cast<double>(j));
					if (!local_smooth && dx > AMR_ratio_*edges_[i] * 1.4*l_reduce)
					{
						smooth = true;
						break;
					}
				}
			}
			if (!smooth)
				continue;
			bool smooth_left = (temp_cells[i + shift].density < temp_cells[i - 1 + shift].density
				* dratio) && (temp_cells[i + shift].density * dratio >
					temp_cells[i - 1 + shift].density)
				&& (temp_cells[i + shift].pressure < temp_cells[i - 1 + shift].pressure
					*pratio) && (temp_cells[i + shift].pressure*pratio >
						temp_cells[i - 1 + shift].pressure);
			bool smooth_right = (temp_cells[i + shift].density <
				temp_cells[i + 1 + shift].density * dratio) && (temp_cells[i + shift].density
					* dratio > temp_cells[i + 1 + shift].density)
				&& (temp_cells[i + shift].pressure < temp_cells[i + 1 + shift].pressure
					*pratio) && (temp_cells[i + shift].pressure*pratio >
						temp_cells[i + 1 + shift].pressure);
			if(smooth_left && smooth_right)
			{
				edge_split.push_back(i);
			}
		}
	}
	if (edge_split.size() > 0)
	{
		size_t Nadd = edge_split.size();
		size_t counter = 0;
		std::vector<double> new_edges(edges_.size() + Nadd);
		std::vector<Primitive> new_cells(N + Nadd);
		std::vector<Extensive> new_extensive(N + Nadd);
		new_edges[0] = edges_[0];
		for (size_t i = 0; i < N; ++i)
		{
			new_edges[i+counter+1] = edges_[i+1];
			new_cells[i + counter] = cells_[i];
			if (edge_split[counter] == i)
			{
				Primitive slope = GetSlope(temp_cells[i + shift - 1], temp_cells[i + shift],
					temp_cells[i + shift + 1], temp_edges[i + shift - 1], temp_edges[i + shift],
					temp_edges[i + shift + 1], temp_edges[i + shift + 2]);
				new_cells[i + counter] = cells_[i] - (0.25*
					(temp_edges[i+shift+1]-temp_edges[i+shift]))*slope;
				new_cells[i + counter].energy = eos_.dp2e(new_cells[i + counter].density,
					new_cells[i + counter].pressure);
				new_cells[i + counter].entropy = eos_.dp2s(new_cells[i + counter].density,
					new_cells[i + counter].pressure);
				new_edges[i + counter+1] = 0.5*(edges_[i] + edges_[i + 1]);
				Extensive old = extensives_[i];
				new_extensive[i+counter] = CalcExtensive(geo_, new_edges, new_cells[i + counter], i + counter);
				++counter;
				new_extensive[i + counter] = old;
				new_extensive[i + counter] -= new_extensive[i + counter-1];
				new_edges[i + counter+1] = edges_[i + 1];
				new_cells[i + counter] = CalcPrimitive(geo_.GetVolume(new_edges, i + counter),
					new_extensive[i + counter], eos_);
				new_cells[i + counter].entropy = eos_.dp2s(new_cells[i + counter].density,
					new_cells[i + counter].pressure);
				new_extensive[i + counter].entropy = new_extensive[i + counter].mass*
					new_cells[i + counter].entropy;
			}
			else
				new_extensive[i] = CalcExtensive(geo_, new_edges, new_cells[i + counter], i + counter);
			if (counter == Nadd && i<(N-1))
			{
				std::copy(edges_.begin() + i + 2, edges_.end(), new_edges.begin() + i + counter + 2);
				std::copy(cells_.begin() + i + 1, cells_.end(), new_cells.begin() + i + counter + 1);
				//std::copy(cells_.begin() + i + 1, cells_.end(), new_extensive.begin() + i + counter + 1);
				break;
			}
		}
		cells_ = new_cells;
		edges_ = new_edges;
		extensives_.resize(cells_.size());
		ReCalcExtensives(cells_);
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
