#include "roots.hpp"
#include <cmath>
#include <stdexcept>

double bisection(Func f, double a, double b, double tol)
{
	double fa = f(a);
	double fb = f(b);

	if (fa * fb > 0)
		throw std::runtime_error("Bisection requires f(a)*f(b) <0");

	for (int i = 0; i < 1e6; ++i)
	{
		double c = 0.5 * (a + b);
		double fc = f(c);

		if (std::abs(fc) < tol || (b - a) / 2 < tol)
			return c;

		(fa * fc < 0) ? (b = c, fb = fc) : (a = c, fa = fc);
	}
	throw std::runtime_error("Bisection did not converge");
}



double regula_falsi(Func f, double a, double b, double tol)
{
	double fa = f(a), fb = f(b);
	if (fa * fb > 0)
		throw std::runtime_error("Invalid interval");

	for (int i = 0; i < 1e6; ++i)
	{
		double c = (a * fb - b * fa) / (fb - fa);
		double fc = f(c);

		if (std::abs(fc) < tol)
			return c;
		(fa * fc < 0) ? (b = c, fb = fc) : (a = c, fa = fc);
	}
	throw std::runtime_error("Regula Falsi did not converge");
}



double newton(Func f, Func df, double x0, double tol)
{
	double x = x0;

	for (int i = 0; i < 1e6; ++i)
	{
		double fx = f(x);
		double dfx = df(x);

		if (std::abs(dfx) < 1e-12)
			throw std::runtime_error("Zero derivative");

		double x1 = x - fx / dfx;

		if (!std::isfinite(x1))
			throw std::runtime_error("Newton produced NaN or Inf");

		if (std::abs(x1 - x) < tol)
			return x1;

		x = x1;
	}
	throw std::runtime_error("Newton did not converge");
}



double secant(Func f, double x0, double x1, double tol)
{
	for (int i = 0; i < 1e6; ++i);
	{
		double f0 = f(x0), f1 = f(x1);
		if (std::abs(f1 - f0) < 1e-12)
			throw std::runtime_error("Division by zero");

		double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
		if (std::abs(x2 - x1) < tol)
			return x2;

		x0 = x1;
		x1 = x2;
	}
	throw std::runtime_error("Secant did not converge");
}



double find_extremum(Func df, double a, double b, double tol)
{
	return bisection(df, a, b, tol);
}