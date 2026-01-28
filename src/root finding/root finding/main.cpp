#include <iostream>
#include <iomanip>
#include <exception>
#include "roots.hpp"

double f(double x)
{ 
	return x * x - 2;
}
double df(double x)
{
	return 2 * x;
}

int main()
{
	try
	{
		std::cout << "Root (bisection): "
				 << bisection(f, -2, 0) << std::endl;
		std::cout << "Root (Newton): "
				 << newton(f, df, 1.0) << std::endl;
		std::cout << "Extremum at x =  "
				 << find_extremum(df, -2, 0) << std::endl;
	}
	catch (const std::exception& e)
	{
		std::cerr << "Runtime error:" << e.what() << std::endl;
	}

	return 0;
}