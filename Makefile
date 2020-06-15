all: parallel_mpi.o parallel_mpi_openmp.o sequential.o

parallel_mpi.o: parallel_mpi.c
	mpicc -std=c99 parallel_mpi.c -o parallel_mpi.o

parallel_mpi_openmp.o: parallel_mpi_openmp.c
	mpicc -std=c99 -fopenmp parallel_mpi_openmp.c -o parallel_mpi_openmp.o

sequential.o: sequential.c
	gcc -std=c99 -fopenmp sequential.c -o sequential.o

clean:
	rm -rf *.o