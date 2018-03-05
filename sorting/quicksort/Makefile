default:
	mpicc qs_mpi.c -o qs_mpi
	mpirun -np 16 ./qs_mpi 10000000 output.txt
	diff output.txt test.txt
