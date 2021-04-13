#include "constants.h"
#include "footingsystem.h"
#include <iostream>
#include <cmath>
#include <iomanip>

double vRdc(double a, const Footingsystem& ref)
{
	double rho1{std::max(0.26 * 0.3 / ref.f_yk * std::pow(ref.f_ck, 2.0/3.0), 0.0013)};
	double p1{(0.18/constants::gcSF) * std::pow(100.0 * rho1, 1.0/3.0)};
	double p2{1.0 + std::sqrt(0.2/ref.d)};
	double p3{0.035 * std::sqrt(p2) * std::pow(ref.f_ck, 1.0/6.0)};
	double p4{pow(ref.f_ck, 1.0/3.0) * 2.0 * ref.d/a};
	double result{std::max(p1, p3) * p2 * p4};
	return 1000.0 * result;
}

double vRd_max(const Footingsystem& ref)
{
	double fcd{0.85 * ref.f_ck / constants::gcSF};
	double result{0.3 * fcd * (1 - ref.f_ck / 250)};
	return 1000.0 * result;
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

double betaFun(double exred, double eyred, double a, double u, double wx, double wy, const Footingsystem& ref)
{
	if (exred > 0 && eyred == 0)
	{
		double ratio{ref.cx / ref.cy};
		if (ratio <= 0.5)
		{
			double k{0.45};
			return 1.0 + k * exred * u * wx;
		}
		else if (ratio == 1.0)
		{
			double k{0.60};
			return 1.0 + k * exred * u * wx;
		}
		else if (ratio == 2.0)
		{
			double k{0.70};
			return 1.0 + k * exred * u * wx;
		}
		else if (ratio >= 3.0)
		{
			double k{0.80};
			return 1.0 + k * exred * u * wx;
		}
	}
	
	else if (exred == 0 && eyred > 0)
	{
		double ratio{ref.cy / ref.cx};
		if (ratio <= 0.5)
		{
			double k{0.45};
			return 1.0 + k * eyred * u * wy;
		}
		else if (ratio == 1.0)
		{
			double k{0.60};
			return 1.0 + k * eyred * u * wy;
		}
		else if (ratio == 2.0)
		{
			double k{0.70};
			return 1.0 + k * eyred * u * wy;
		}
		else if (ratio >= 3.0)
		{
			double k{0.80};
			return 1.0 + k * eyred * u * wy;
		}
	}
	
	else if (exred > 0 && eyred > 0)
	{
		double p1{std::pow(exred / (ref.cy + 2.0 * a), 2)};
		double p2{std::pow(eyred / (ref.cx + 2.0 * a), 2)};
		double p3{std::sqrt(p1 + p2)};
		return 1.0 + 1.8 * p3;
	}
}

double limitFinder(const Footingsystem& ref)
{
	double pos_x{ref.bx / 2.0 - ref.ax - ref.cx / 2.0};
	double neg_x{ref.bx / 2.0 + ref.ax - ref.cx / 2.0};
	double pos_y{ref.by / 2.0 - ref.ay - ref.cy / 2.0};
	double neg_y{ref.bx / 2.0 + ref.ax - ref.cx / 2.0};
	double min_x{std::min(pos_x, neg_x)};
	double min_y{std::min(pos_y, neg_y)};
	return std::min(min_x, min_y);
}

void printer(const Footingsystem& ref)
{
	double a{0.0};
	const double* ptrIncr{&constants::step};
	double limit_equation{2.0 * ref.d};
	double limit_geometry{limitFinder(ref)};
	
	while(a <= limit_equation)
	{
		
		std::cout << "Distance from column: [a = " << a << "] --> ";
		
		if (a <= limit_geometry)
			std::cout << "(OK)" << '\n';
		else 
			std::cout << "(NOT OK)" << '\n';
	
		double temp_medxred{MEd_x_red(e_x(ref), e_y(ref), area(a, ref), weight(ref), I_x(a, ref), ref)};
		double temp_medyred{MEd_y_red(e_x(ref), e_y(ref), area(a, ref), weight(ref), I_y(a, ref), ref)};
		double temp_vedred{VEd_red(e_x(ref), e_y(ref), area(a, ref), weight(ref), ref)};
		double temp_exred{e_x_red(temp_medxred, temp_vedred)};
		double temp_eyred{e_y_red(temp_medyred, temp_vedred)};
		double temp_beta{betaFun(temp_exred, temp_eyred, a, u(a, ref), W_x(a, ref), W_y(a, ref), ref)};
		double temp_maxved{max_vEd(temp_vedred, u(a, ref), temp_beta, ref)};
		double temp_vrdc{vRdc(a, ref)};
		
		std::cout << "u(a) 		     = " << u(a, ref) << '\n';
		std::cout << "A'(a)		     = " << area(a, ref) << '\n';
		std::cout << "Ved,red(a)  	     = " << temp_vedred << '\n';
		std::cout << "Ix'(a)		     = " << I_x(a, ref) << '\n';
		std::cout << "Iy'(a)		     = " << I_y(a, ref) << '\n';
		//std::cout << "Wx 		     = " << W_x(a, ref) << '\n';
		//std::cout << "Wy 		     = " << W_y(a, ref) << '\n';
		std::cout << "MEd,xred(a)          = " << temp_medxred << '\n';
		std::cout << "MEd,yred(a)          = " << temp_medyred << '\n';
		std::cout << "ex,red(a)   	     = " << temp_exred << '\n';
		std::cout << "ey,red(a)   	     = " << temp_eyred << '\n';
		std::cout << "beta(a)     	     = " << temp_beta << '\n';
		std::cout << "vRd,c(a)    	     = " << temp_vrdc << '\n';
		std::cout << "maxvEd(a)   	     = " << temp_maxved << '\n';
		std::cout << "vRd,c(a) / maxvEd(a) = " << temp_vrdc / temp_maxved << '\n';
		//std::cout << "vRd,max 	     = " << vRd_max(ref) << '\n';

		std::cout << "----------------------------------------------------" << '\n';
		a += *ptrIncr;
	}

	ptrIncr = nullptr;
}