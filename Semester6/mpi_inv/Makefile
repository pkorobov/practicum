all: mpi_inv
mpi_inv: main.o mpi_algo.o mpi_etc.o mpi_residual.o
	mpicxx main.o mpi_algo.o mpi_etc.o mpi_residual.o -o mpi_inv -O3

main.o: main.cpp
	mpicxx -c main.cpp -O3

mpi_algo.o: mpi_algo.cpp
	mpicxx -c mpi_algo.cpp -O3

mpi_etc.o: mpi_etc.cpp
	mpicxx -c mpi_etc.cpp -O3

mpi_residual.o: mpi_residual.cpp
	mpicxx -c mpi_residual.cpp -O3

clean:
	rm -rf *.o mpi_inv 
