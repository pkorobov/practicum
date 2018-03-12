all: mpi_inv
mpi_inv: main.o mpi_algo.o mpi_etc.o
	g++ main.o mpi_algo.o mpi_etc.o -o hello

main.o: main.cpp
	g++ -c main.cpp

mpi_algo.o: mpi_algo.cpp
	g++ -c mpi_algo.cpp

mpi_etc.o: mpi_etc.cpp
	g++ -c mpi_etc.cpp

clean:
	rm -rf *.o mpi_inv 
