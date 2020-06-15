all: parallel_mpi.o parallel_mpi_openmp.o sequential.o

parallel_mpi.o: parallel_mpi.c
	mpicc parallel_mpi.c -o parallel_mpi.o

parallel_mpi_openmp.o: parallel_mpi_openmp.c
	mpicc -fopenmp parallel_mpi_openmp.c -o parallel_mpi_openmp.o

sequential.o: sequential.c
	gcc -fopenmp -std=c99 sequential.c -o sequential.o

clean:
	rm -rf *.o
