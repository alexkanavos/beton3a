#include "footingsystem.h"
#include "calculations.h"

int main()
{
	Footingsystem sample{2.3, 2.3, 0.6, 0.55, 0.6, 0.4, 0.0, 0.0, 1.0, 25000.0, 500000.0, 550.0, 180.0, 1200.0};
	
	printer(sample);
	
	return 0;
}