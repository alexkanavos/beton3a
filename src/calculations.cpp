#include "constants.h"
#include "footingsystem.h"
#include <iostream>
#include <cmath>

void initializer(Footingsystem& ref)
{
	std::cout << "Insert required input values. \n";

	std::cout << "bx: ";
    std::cin >> ref.bx;
    
    std::cout << "by: ";
    std::cin >> ref.by;
    
    std::cout << "h: ";
    std::cin >> ref.h;

    std::cout << "d: ";
    std::cin >> ref.d;

    std::cout << "cx: ";
    std::cin >> ref.cx;
    
    std::cout << "cy: ";
    std::cin >> ref.cy;
    
    std::cout << "ax: ";
    std::cin >> ref.ax;
    
    std::cout << "ay: ";
    std::cin >> ref.ay;

    std::cout << "t: ";
    std::cin >> ref.t;
    
    std::cout << "Mx: ";
    std::cin >> ref.Mx;
    
    std::cout << "My: ";
    std::cin >> ref.My;

    std::cout << "Ntot: ";
    std::cin >> ref.Ntot;
}

double vRdc(double a, const Footingsystem& ref)
{
	double p1{(0.18/1.5) * std::pow(100.0 * constants::rho1, 1.0/3.0)};
	double p2{1.0 + std::sqrt(0.2/ref.d)};
	double p3{0.035 * std::sqrt(p2) * std::pow(constants::f_ck, 1.0/6.0)};
	double p4{pow(constants::f_ck, 1.0/3.0) * 2 * ref.d/a};
	double result{std::max(p1, p3) * p2 * p4};
	return result;
}

double area(double a, const Footingsystem& ref)
{
	double result{constants::pi * std::pow(a, 2) + ref.cx * ref.cy + 2.0 * a * (ref.cx + ref.cy)};
	return result;
}

double weight(const Footingsystem& ref)
{
	double result{ref.bx * ref.by * ref.h * constants::g_con + constants::g_soil * (ref.bx * ref.by - ref.cx * ref.cy) * (ref.t - ref.h)};
	return result;
}

double e_x(const Footingsystem& ref)
{
	return ref.Mx / ref.Ntot;
}

double e_y(const Footingsystem& ref)
{
	return ref.My / ref.Ntot;
}

double VEd_red(double ex, double ey, double area, double weight, const Footingsystem& ref)
{
	double p1{1.0 - area/(ref.bx * ref.by)};
	double p2{ref.Ntot - weight};
	double p3{12.0 * area/(ref.bx * ref.by)};
	double p4{(ref.ax * ex / std::pow(ref.bx, 2) + ref.ay * ey / std::pow(ref.by, 2)) * ref.Ntot};
	double result{p1 * p2 - p3 * p4};
	return result;
}

double u(double a, const Footingsystem& ref)
{
	double result{2.0 * (constants::pi * a + ref.cx + ref.cy)};
	return result;
}

double max_vEd(double vedred, double u, double beta, const Footingsystem& ref)
{
	double result{beta * vedred / (u * ref.d)};
	return result;
}

double I_x(double a, const Footingsystem& ref)
{
	double p1{(1.0/4.0) * constants::pi * std::pow(a, 2) * (std::pow(ref.cx, 2) + std::pow(a, 2))};
	double p2{(1.0/6.0) * a * std::pow(ref.cx, 3)};
	double p3{(1.0/12.0) * ref.cy * std::pow(2 * a + ref.cx, 3)};
	double result{p1 + p2 + p3};
	return result;
}

double I_y(double a, const Footingsystem& ref)
{
	double p1{(1.0/4.0) * constants::pi * std::pow(a, 2) * (std::pow(ref.cy, 2) + std::pow(a, 2))};
	double p2{(1.0/6.0) * a * std::pow(ref.cy, 3)};
	double p3{(1.0/12.0) * ref.cx * std::pow(2.0 * a + ref.cy, 3)};
	double result{p1 + p2 + p3};
	return result;
}

double W_x(double a, const Footingsystem& ref)
{
	double result{(1.0/2.0) * std::pow(ref.cx, 2) + ref.cx * ref.cy + 2.0 * a * ref.cy + constants::pi * a * ref.cx + 4.0 * std::pow(a, 2)};
	return result;
}

double W_y(double a, const Footingsystem& ref)
{
	double result{(1.0/2.0) * std::pow(ref.cy, 2) + ref.cx * ref.cy + 2.0 * a * ref.cx + constants::pi * a * ref.cy + 4.0 * std::pow(a, 2)};
	return result;
}

double MEd_x_red(double ex, double ey, double area, double weight, double giotax, const Footingsystem& ref)
{
	double p1{ex * ref.Ntot};
	double p2{area * ref.ax * (ref.Ntot - weight) / (ref.bx * ref.by)};
	double p3{12.0 * ref.Ntot / (ref.bx * ref.by)};
	double p4{ex * (giotax + area * std::pow(ref.ax, 2)) / std::pow(ref.bx, 2)};
	double p5{ey * ref.ax * ref.ay * area / std::pow(ref.by, 2)};
	double result{p1 - p2 - p3 * (p4 + p5)};
	return result;
}

double MEd_y_red(double ex, double ey, double area, double weight, double giotay, const Footingsystem& ref)
{
	double p1{ey * ref.Ntot};
	double p2{area * ref.ay * (ref.Ntot - weight) / (ref.bx * ref.by)};
	double p3{12.0 * ref.Ntot / (ref.bx * ref.by)};
	double p4{ey * (giotay + area * std::pow(ref.ay, 2)) / std::pow(ref.by, 2)};
	double p5{ex * ref.ax * ref.ay * area / std::pow(ref.bx, 2)};
	double result{p1 - p2 - p3 * (p4 + p5)};
	return result;
}

double e_x_red(double medxred, double vedred)
{
	return medxred / vedred;
}

double e_y_red(double medyred, double vedred)
{
	return medyred / vedred;
}

double betaFun(double exred, double eyred, double u, double wx, double wy, double a, const Footingsystem& ref)
{
	double result{};
	if (exred > 0 && eyred == 0)
	{
		double logos{ref.cx / ref.cy};
		double k{};
		double* ptr{&k};
		if (logos <= 0.5)
		{
			*ptr = 0.45;
		}
		else if (logos = 1.0)
		{
			*ptr = 0.60;
		}
		else if (logos = 2.0)
		{
			*ptr = 0.70;
		}
		else if (logos >= 3.0)
		{
			*ptr = 0.80;
		}
		
 		result = 1.0 + k * exred * u / wx;
 	}
 	else if (exred == 0 && eyred > 0)
 	{
 		double logos{ref.cy / ref.cx};
 		double k{};
		double* ptr{&k};
		if (logos <= 0.5)
		{
			*ptr = 0.45;
		}
		else if (logos = 1.0)
		{
			*ptr = 0.60;
		}
		else if (logos = 2.0)
		{
			*ptr = 0.70;
		}
		else if (logos >= 3.0)
		{
			*ptr = 0.80;
		}
		
 		result = 1.0 + k * eyred * u / wy;
 	}
    else if (exred > 0 && eyred > 0)
    {
    	double p1{std::pow(exred/(ref.cy + 2.0 * a), 2)};
    	double p2{std::pow(eyred/(ref.cx + 2.0 * a), 2)};
    	double p3{std::sqrt(p1 + p2)};
    	return 1.0 + 1.8 * p3;
    }
    return result;
}