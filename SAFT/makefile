saft: saft.o bscan.o plate.o img.o
	g++ -o saft saft.o bscan.o plate.o img.o

saft.o: saft.cpp saft.h
	g++ -c saft.cpp
bscan.o: bscan.cpp saft.h
	g++ -c bscan.cpp
plate.o: plate.cpp saft.h
	g++ -c plate.cpp
img.o: img.cpp saft.h
	g++ -c img.cpp
