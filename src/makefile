output: main.o calculations.o
	g++ -std=c++17 main.o calculations.o -o output

main.o: main.cpp calculations.h footingsystem.h
	g++ -std=c++17 -c main.cpp

calculations.o: calculations.cpp calculations.h constants.h footingsystem.h
	g++ -std=c++17 -c calculations.cpp

clean:
	rm *.o output