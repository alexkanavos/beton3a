#include "footingsystem.h"
#include "calculations.h"
#include <iostream>

int main()
{
	Footingsystem sample{};
	
	initializer(sample);
	
	double test{vRdc(0.9, sample)};

	std::cout << test << '\n';
}