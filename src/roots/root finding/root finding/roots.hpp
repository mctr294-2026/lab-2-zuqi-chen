#ifndef ROOTS_HPP
#define ROOTS_HPP

using Func = double (*)(double);

double bisection(Func f, double a, double b, double tol = 1e-6);

double regula_falsi(Func f, double a, double b, double tol = 1e-6);

double newton(Func f, Func df, double x0, double tol = 1e-6);

double secant(Func f, double x0, double x1, double tol = 1e-6);

double find_extremum(Func df, double a, double b, double tol = 1e-6);

#endif