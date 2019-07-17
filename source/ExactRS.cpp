#include "ExactRS.hpp"
#include "universal_error.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>
#if defined(_MSC_VER)
/* Microsoft C/C++-compatible compiler */
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

namespace
{
	double fastsqrt(double x)
	{
		if (x<std::numeric_limits<float>::min() || x > std::numeric_limits<float>::max())
			return std::sqrt(x);
		double res = static_cast<double>(_mm_cvtss_f32(_mm_rsqrt_ss(_mm_set_ss(static_cast<float>(x)))));
		return x*res*(1.5 - 0.5*res*res*x);
	}

	double CalcFrarefraction(Primitive const& cell, double p,double gamma)
	{
		double cs = fastsqrt(gamma*cell.pressure/ cell.density);
		return 2 * cs*(std::pow(p / cell.pressure, (gamma - 1) / (2 * gamma)) - 1) / (gamma - 1);
	}

	double CalcFshock(Primitive const& cell, double p,double gamma)
	{
		double A = 2 / ((gamma + 1) *cell.density);
		double B = (gamma - 1)*cell.pressure / (gamma + 1);
		return (p - cell.pressure)*fastsqrt(A / (p + B));
	}

	double dCalcFshock(Primitive const& cell, double p, double gamma)
	{
		double A = 2 / ((gamma + 1)*cell.density);
		double B = (gamma - 1)*cell.pressure / (gamma + 1);
		return fastsqrt(A / (B + p))*(1 - (p - cell.pressure) / (2 * (B + p)));
	}

	double dCalcFrarefraction(Primitive const& cell, double p, double gamma)
	{	
		double cs = fastsqrt(gamma*cell.pressure/ cell.density);
		return std::pow(cell.pressure / p, (gamma + 1) / (2 * gamma)) / (cell.density*cs);
	}

	RSsolution GetFirstGuess(Primitive const & left, Primitive const & right,double gamma)
	{
		double Pmax = std::max(left.pressure, right.pressure);
		double Pmin = std::min(left.pressure, right.pressure);
		double Q = Pmax / Pmin;
		double csl = fastsqrt(gamma*left.pressure / left.density);
		double csr = fastsqrt(gamma*right.pressure / right.density);
		RSsolution res;
		// First try Ppvrs
		double pvrs = 0.5*(Pmax + Pmin) + 0.125*(left.velocity - right.velocity)*(left.density + right.density)*
			(csl + csr);
		if (Q < 2 && pvrs >= Pmin && pvrs <= Pmax)
		{
			res.pressure = pvrs;
			return res;
		}
			
		if (pvrs <= Pmin) // Use Two rarefactions
		{
			double z = (gamma - 1) / (2 * gamma);
			res.pressure = std::pow((csl + csr - (gamma - 1)*(right.velocity - left.velocity)*0.5) / (csl*
				std::pow(left.pressure, -z)	+ csr*std::pow(right.pressure, -z)), 1 / z);
		}
		else // Two shocks
		{
			double Al = 2 / ((gamma + 1)*left.density);
			double Ar = 2 / ((gamma + 1)*right.density);
			double Bl = (gamma - 1)*left.pressure / (gamma + 1);
			double Br = (gamma - 1)*right.pressure / (gamma + 1);
			res.pressure = std::max(0.0, res.pressure);
			double gl = fastsqrt(Al / (res.pressure + Bl));
			double gr = fastsqrt(Ar / (res.pressure + Br));
			res.pressure = (gl*left.pressure + gr*right.pressure + left.velocity - right.velocity) / (gl + gr);
			if (res.pressure < Pmin)
				res.pressure = pvrs;
		}
		return res;
	}

	double GetValue(Primitive const & left, Primitive const & right,double p,double gamma)
	{
		double res = right.velocity-left.velocity;
		if (p > left.pressure)
			res += CalcFshock(left, p, gamma);
		else
			res += CalcFrarefraction(left, p,gamma);
		if (p > right.pressure)
			res += CalcFshock(right, p, gamma);
		else
			res += CalcFrarefraction(right, p, gamma);
		return res;
	}

	double GetdValue(Primitive const & left, Primitive const & right, double p, double gamma)
	{
		double res = 0;
		if (p > left.pressure)
			res += dCalcFshock(left, p, gamma);
		else
			res += dCalcFrarefraction(left, p, gamma);
		if (p > right.pressure)
			res += dCalcFshock(right, p, gamma);
		else
			res += dCalcFrarefraction(right, p, gamma);
		return res;
	}

	double Bisection(Primitive const& left, Primitive const& right,double gamma,
		double max_scale,double minp,double guess)
	{
		double a = guess*0.2;
		double b = guess * 5;
		double c = guess;
		double dp = 0;
		double valuec = GetValue(left, right, b , gamma);
		double valuea = GetValue(left, right, a, gamma);
		if (valuea*valuec > 0)
		{
			a = guess * 0.0001;
			b = guess * 5;
			c = guess;
			dp = 0;
			valuec = GetValue(left, right, b, gamma);
			valuea = GetValue(left, right, a, gamma);
			if (valuea*valuec > 0)
				throw UniversalError("Same sign in RS");
		}
			
		int counter = 0;
		while (((b-a) > 1e-10*(minp + c)) || (std::abs(valuec) > 1e-7*max_scale))
		{
			c = (a + b)*0.5;
			valuec = GetValue(left, right, c, gamma);
			if (valuec*valuea > 0)
			{
				a = c;
				valuea = valuec;
			}
			else
				b = c;
			++counter;
			if (counter > 300)
			{
				UniversalError eo("Too many iterations in RS");
				eo.AddEntry("Left density", left.density);
				eo.AddEntry("Left pressure", left.pressure);
				eo.AddEntry("Left velocity", left.velocity);
				eo.AddEntry("Right density", right.density);
				eo.AddEntry("Right pressure", right.pressure);
				eo.AddEntry("Right velocity", right.velocity);
				throw eo;
			}
		}
		return c;
	}
}

ExactRS::ExactRS(double gamma):gamma_(gamma)
{}

ExactRS::~ExactRS()
{}

RSsolution ExactRS::Solve(Primitive const & left, Primitive const & right)const
{
	const double eps = 1e-7;
	// Is there a vaccum?
	double dv = right.velocity - left.velocity;
	double soundspeeds = 2 * (fastsqrt(gamma_*left.pressure / left.density) + fastsqrt(gamma_*right.pressure / right.density)) / (gamma_ - 1);
	if (dv >= soundspeeds)
	{
		RSsolution res;
		res.pressure = 0;
		res.velocity = 0.5*(right.velocity+left.velocity);
		return res;
	}
	RSsolution res = GetFirstGuess(left, right,gamma_);
	res.velocity = 0;
	if (res.pressure < 0)
	{
		res.pressure = 0;
		return res;
	}
	double value = GetValue(left, right, res.pressure, gamma_);
	double Cs = std::max(fastsqrt(gamma_*left.pressure / left.density), fastsqrt(gamma_*right.pressure / right.density));
	double max_scale = std::max(Cs, std::max(std::abs(left.velocity), std::abs(right.velocity)));
	double dp = 0;
	double minp = std::min(left.pressure, right.pressure);
	double p = res.pressure;
	size_t counter = 0;
	while ((abs(dp) > eps*(minp + res.pressure)) || (std::abs(value) > eps*max_scale))
	{
		p = res.pressure;
		dp = value / GetdValue(left, right, res.pressure, gamma_);
		res.pressure -= std::max(std::min(dp,res.pressure*0.5),-0.5*res.pressure);
		value = GetValue(left, right, res.pressure, gamma_);
		++counter;
		if (counter > 30)
		{
			res.pressure = Bisection(left, right, gamma_, max_scale, minp,res.pressure);
			break;
		}
	} 
	double fr = (res.pressure > right.pressure) ? CalcFshock(right, res.pressure, gamma_) : CalcFrarefraction(
		right, res.pressure, gamma_);
	double fl = (res.pressure > left.pressure) ? CalcFshock(left, res.pressure, gamma_) : CalcFrarefraction(
		left, res.pressure, gamma_);
	res.velocity = 0.5*(left.velocity + right.velocity) + 0.5*(fr-fl);
	return res;
}
