#ifndef CALCULATIONS_H
#define CALCULATIONS_H

void initializer(Footingsystem& ref);

double vRdc(double a, const Footingsystem& ref);

double area(double a, const Footingsystem& ref);

double weight(const Footingsystem& ref);

double e_x(const Footingsystem& ref);

double e_y(const Footingsystem& ref);

double VEd_red(double ex, double ey, double area, double weight, const Footingsystem& ref);

double u(double a, const Footingsystem& ref);

double max_vEd(double vedred, double u, double beta, const Footingsystem& ref);

double I_x(double a, const Footingsystem& ref);

double I_y(double a, const Footingsystem& ref);

double W_x(double a, const Footingsystem& ref);

double W_y(double a, const Footingsystem& ref);

double MEd_x_red(double ex, double ey, double area, double weight, double giotax, const Footingsystem& ref);

double MEd_y_red(double ex, double ey, double area, double weight, double giotay, const Footingsystem& ref);

double e_x_red(double medxred, double vedred);

double e_y_red(double medyred, double vedred);

double betaFun(double exred, double eyred, double u, double wx, double wy, double a, const Footingsystem& ref);

#endif