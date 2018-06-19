#define _USE_MATH_DEFINES
#include "hdsim.hpp"
#include "hdf_util.hpp"
#include "MinMod.hpp"
#include <iostream>
#include <fstream>
#include <cassert>
#include <boost/math/tools/roots.hpp>
#include <float.h>
#include <Windows.h>

unsigned int oldd = 0;
unsigned int fp_control_state = _controlfp_s(&oldd,_EM_INEXACT, _MCW_EM);

//#define restart

namespace
{
	string int2str(int n)
	{
		stringstream ss;
		ss << n;
		return ss.str();
	}

	template <typename T>
	T LinearInterpolation(const vector<T> &x, const vector<T> &y, T xi)
	{
		typename vector<T>::const_iterator it = upper_bound(x.begin(), x.end(), xi);
		assert(it != x.end() && it != x.begin() &&
			"X out of range in Linear Interpolation");
		if (*it == xi)
			return y[static_cast<size_t>(it - x.begin())];

		// Are we near the edge?
		if (it == x.end() - 1)
			return y[static_cast<size_t>(it - x.begin())] + (xi - *it)*
			(y[static_cast<size_t>(it - 1 - x.begin())] -
				y[static_cast<size_t>(it - x.begin())]) / (*(it - 1) - *it);
		else
			return y[static_cast<size_t>(it - x.begin())] + (xi - *it)*
			(y[static_cast<size_t>(it + 1 - x.begin())] -
				y[static_cast<size_t>(it - x.begin())]) / (*(it + 1) - *it);
	}

	bool missing_file_data(string const& fname)
	{
		std::cout << "Could not find file " << fname << std::endl;
		return false;
	}

	vector<double> read_vector(string const& fname)
	{
		double buf = 0;
		vector<double> res;
		ifstream f(fname.c_str());
		assert(f || missing_file_data(fname));
		while (f >> buf)
			res.push_back(buf);
		f.close();
		return res;
	}

	vector<Primitive> calc_init(vector<double> const& edges, double M, double R, double Nemden,double stargamma)
	{
		// Get lane emden result, G=1 R=Rsun M=Msun
		vector<double> xsi, theta, dtheta;
		if (stargamma < 1.6)
		{
			xsi = read_vector("c:/xsin3.txt");
			theta = read_vector("c:/thetan3.txt");
			dtheta = read_vector("c:/dthetan3.txt");
		}
		else
		{
			xsi = read_vector("c:/xsi.txt");
			theta = read_vector("c:/theta.txt");
			dtheta = read_vector("c:/dtheta.txt");
		}
		size_t n = xsi.size();
		double rhoc = -M*xsi[n - 1] / (4 * M_PI*R*R*R*dtheta[n - 1]);
		double K = M*M*pow(4 * M_PI, 1.0 / Nemden)*pow(-M*xsi[n - 1] / (dtheta[n - 1] * R*R*R), -(Nemden + 1.0) / Nemden) /
			((dtheta[n - 1] * dtheta[n - 1])*(1 + Nemden)*pow(R, 4));
		vector<Primitive> res(edges.size() - 1);
		for (size_t i = 0; i < res.size(); ++i)
		{
			double y = 0.5 * xsi[n - 1] * (edges[i + 1] + edges[i]) / R;
			double Theta = 0;
			if (y > xsi[n - 1])
			{
				throw "bad size";
			}
			else
				Theta = LinearInterpolation(xsi, theta, y);
			res[i].density = rhoc*pow(Theta, Nemden);
			res[i].pressure = K*pow(res[i].density, (Nemden + 1) / Nemden);
			res[i].velocity = 0;
		}
		return res;
	}

	vector<double> getedges(size_t N, double L)
	{
		double factor = 0.9999;
		double R = L*factor;
		size_t Nouter = (1 - factor) * 5 * N;
		size_t Nmid = N - Nouter;
		double dx = R / Nmid;
		vector<double> res(N + 1);
		for (size_t i = 0; i <= Nmid; ++i)
			res[i] = i*dx;
		dx = (L - res[Nmid])/(Nouter+1);
		for (size_t i = Nmid+1; i <= N; ++i)
			res[i] = res[i - 1] + dx;
		return res;
	}

	class ParaboleAnomaly
	{
	private:
		double t_, Rp_, Mbh_;
	public:
		ParaboleAnomaly(double t, double Rp, double Mbh) :t_(t), Rp_(Rp), Mbh_(Mbh) {}
		double operator()(double f)
		{
			return t_ - sqrt(2 * pow(Rp_, 3) / Mbh_)*tan(f / 2)*(3 + pow(tan(f / 2), 2)) / 3;
		}
	};

	struct TerminationCondition
	{
		bool operator() (double min, double max)
		{
			return abs(min - max) <= 0.000001;
		}
	};

	double GetTrueAnomaly(double t, double Rp, double Mbh)
	{

		ParaboleAnomaly Eanom(t, Rp, Mbh);
		pair<double, double> res = boost::math::tools::bisect(Eanom, -M_PI, M_PI,
			TerminationCondition());
		return res.first;
	}

	class Gravity : public SourceTerm
	{
	private:
		vector<double> acc_;
		mutable vector<double> total_acc_;
		double rhoc_;
		double alpha_;
		double Mbh_;
		double Rp_;
		mutable double f_;
		bool selfgravity_;
	public:
		Gravity(double M, double R, double Mbh, double Rp,bool selfgravity,double stargamma, vector<double> const& edges) :
			acc_(vector<double>()), total_acc_(vector<double>()),rhoc_(0),alpha_(0), Mbh_(Mbh), Rp_(Rp), f_(0),selfgravity_(selfgravity)
		{
			vector<double> xi, dtheta;
			if (stargamma > 1.6)
			{
				xi = read_vector("c:/xsi.txt");
				dtheta = read_vector("c:/dtheta.txt");
			}
			else
			{
				xi = read_vector("c:/xsin3.txt");
				dtheta = read_vector("c:/dthetan3.txt");
			}
			size_t n = xi.size();
			rhoc_ = -M*xi[n - 1] / (4 * M_PI*R*R*R*dtheta[n - 1]);
			alpha_ = R / xi[n - 1];

			size_t N = edges.size()-1;
			acc_.resize(N);
			for (size_t i = 0; i < N; ++i)
			{
				double x = 0.5*(edges[i + 1] + edges[i]);
				double xsi = x / alpha_;
				if (xsi > xi.back())
					xsi = xi.back() - 0.000001;
				double Mtemp = -rhoc_*alpha_*alpha_*alpha_ * 4 * M_PI*LinearInterpolation(xi, dtheta, xsi)*xsi*xsi;
				acc_[i] = -Mtemp / (x*x);
			}
			total_acc_ = acc_;
		}

		void CalcForce(vector<double> const& edges, vector<Primitive> const& cells, double time,
			vector<Extensive> & extensives, double dt)const
		{
			size_t N = cells.size();
			f_ = GetTrueAnomaly(time, Rp_, Mbh_);
			for (size_t i = 0; i < N; ++i)
			{
				double acc = acc_[i];
				if (!selfgravity_)
					acc = 0;
				// Do Mbh
				double R = 2 * Rp_ / (1 + cos(f_));
				double x = 0.5*(edges[i + 1] + edges[i]);
				acc -= 0*Mbh_*x*pow(sqrt(R*R + x*x), -3);
				total_acc_[i] = acc;
				extensives[i].momentum += extensives[i].mass*acc*dt;
				extensives[i].energy += extensives[i].mass*acc*dt*cells[i].velocity;
			}
		}

		double GetInverseTimeStep(vector<double> const& edges)const
		{
			double res = 0;
			size_t Nloop = total_acc_.size();
			for (size_t i = 0; i < Nloop; ++i)
				res = max(res, sqrt(abs(total_acc_[i])/ (edges[i + 1] - edges[i])));
			return res;
		}
	};

	double read_number(string const& fname)
	{
		double buf = 0;
		ifstream f(fname.c_str());
		assert(f || missing_file_data(fname));
		f >> buf;
		f.close();
		return buf;
	}

	void ReadInput(double &beta,double &star_gamma,double &gamma,bool &selfgravity)
	{
		beta = read_number("c:/beta.txt");
		star_gamma = read_number("c:/star_gamma.txt");
		gamma = read_number("c:/gamma.txt");
		double sg = read_number("c:/selfgravity.txt");
		selfgravity = sg > 0.5;
	}

	void DragForce(hdsim &sim)
	{
		vector<Extensive> &extensive = sim.GetExtensives();
		vector<Primitive> &cells = sim.GetCells();
		size_t N = cells.size();
		for (size_t i = 0; i < N; ++i)
		{
			cells[i].velocity *= 0.99999;
			double dp = 0.01*extensive[i].momentum;
			double p = extensive[i].momentum;
			extensive[i].momentum *= 0.99999;
			extensive[i].energy =extensive[i].et + 0.5*extensive[i].momentum*extensive[i].momentum / extensive[i].mass;
		}
	}
}

int main(void)
{
	// Units G=1 M=solar R=solar t=1.592657944577715e+03
	bool selfgravity = false;
	double star_gamma = 4. / 3.;
	double gamma = 4./3.;
	double beta = 8;
	ReadInput(beta, star_gamma, gamma, selfgravity);
	double cfl = 0.2;
	if (beta < 7)
		cfl = 0.5;
	else
		if (beta < 15)
			cfl = 0.3;
	ExactRS rs(gamma);
	IdealGas eos(gamma);
	RigidWall bl;
	//ConstantPrimitive br(Primitive(1e-25, 1e-26, 0, eos.dp2s(1e-25, 1e-26)));
	//Ratchet br;
	RigidWall br;
	SeveralBoundary boundary(bl, br);
	MinMod interp(boundary);
	double R = 1;
	double M = 1;
	double Mbh = 1e6;
	double Rt = R*pow(Mbh / M, 1.0 / 3.0);
	double Rp = Rt / beta;
#ifdef restart
	Snapshot snap = read_hdf5_snapshot("c:/TidalData/gamma43/gas43/sg/beta129/tide_114.h5");
	vector<double> edges = snap.edges;
	vector<Primitive> cells = snap.cells;
	int counter = 115;
#else
	int counter = 0;
	size_t Np = 512*4;
	if (beta > 30)
		Np = static_cast<size_t>(Np * beta / 30);
	double Nemden = 1 / (star_gamma - 1);
	double fstart = -acos(2 / beta - 1);
	double tstart = sqrt(2 * pow(Rt / beta, 3) / Mbh)*tan(fstart / 2)*(3 + pow(tan(fstart / 2), 2)) / 3;
	vector<double> edges = getedges(Np, R);
	vector<Primitive> cells = calc_init(edges, M, R, Nemden, star_gamma);
#endif

	Gravity source(M, R, Mbh, Rp,selfgravity,star_gamma,edges);
	hdsim sim(cfl, cells, edges, interp, eos, rs,source);
#ifdef restart
	sim.SetTime(snap.time);
	sim.SetCycle(snap.cycle);
#else
	sim.SetTime(tstart);
#endif

	double dt = 0.05;
	double initd =  cells[0].density;
	double maxd = cells[0].density;
	double last = sim.GetTime();
	double mind = maxd;
	initd *= selfgravity ? 1 : 0.8;
	while (sim.GetCells()[0].density> max(0.3*initd,0.2*maxd) && sim.GetTime()<5.52)
	{
		if (sim.GetCycle() % 500 == 0)
			write_snapshot_to_hdf5(sim, "c:/sim_data/temp.h5");
		if (sim.GetCycle() % 100 == 0)
			cout << "Time = " << sim.GetTime() << " Cycle = " << sim.GetCycle() << endl;
		if (sim.GetTime() - last > dt || sim.GetCycle() == 0 || sim.GetCells()[0].density>1.02*maxd || 
			sim.GetCells()[0].density*1.02<mind || last*sim.GetTime()<0)
		{
			string dirloc = "c:\\TidalData";
			CreateDirectory(dirloc.c_str(), NULL);
			if (star_gamma > 1.6)
				dirloc += "\\gamma53";
			else
				dirloc += "\\gamma43";
			CreateDirectory(dirloc.c_str(), NULL);
			if (gamma>1.7)
				dirloc += "\\gas2";
			else
				if (gamma > 1.6)
					dirloc += "\\gas53";
				else
					if(gamma>1.4)
						dirloc += "\\gas15";
					else
						dirloc += "\\gas43";
			CreateDirectory(dirloc.c_str(), NULL);
			dirloc += (selfgravity) ? "\\sg" : "\\nosg";
			CreateDirectory(dirloc.c_str(), NULL);
			dirloc += "\\beta";
			dirloc += int2str(int(beta));
			CreateDirectory(dirloc.c_str(), NULL);

			if(last*sim.GetTime()<0)
				write_snapshot_to_hdf5(sim, dirloc + "\\special.h5");
			else
			{
				write_snapshot_to_hdf5(sim, dirloc + "\\tide_" + int2str(counter) + ".h5");
				++counter;
			}
			last = sim.GetTime();
			maxd = max(maxd, sim.GetCells()[0].density);
			mind = sim.GetCells()[0].density;
		}
		sim.TimeAdvance2();
		DragForce(sim);
	}
//	write_snapshot_to_hdf5(sim, "c:/sim_data/snap1d.h5");
	return 0;
}