output: main.o calculations.o
	g++ -std=c++17 -Wall main.o calculations.o -o output

main.o: main.cpp calculations.h footingsystem.h
	g++ -std=c++17 -Wall -c main.cpp

calculations.o: calculations.cpp calculations.h constants.h footingsystem.h
	g++ -std=c++17 -Wall -c calculations.cpp

clean:
	rm *.o output