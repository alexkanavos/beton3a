#include "footingsystem.h"
#include "calculations.h"

int main()
{
	//Struct initialization
	Footingsystem sample{2.3, 2.3, 0.6, 0.55, 0.6, 0.4, 0.0, 0.0, 1.0, 25.0, 500.0, 550.0, 180.0, 1200.0};

	//Call printer function to get the output
	printer(sample);

	//exit point
	return 0;
}