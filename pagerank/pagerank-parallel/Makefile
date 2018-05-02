pagerank: pagerank-parallel.c pr_graph.c
	mpicc pagerank-parallel.c -o pagerank-parallel			
	mpirun -np 1 ./pagerank-parallel A.graph A.pr.txt	
	mpirun -np 1 ./pagerank-parallel B.graph B.pr.txt	
	mpirun -np 1 ./pagerank-parallel live-journal.graph live-journal.pr.txt
	mpirun -np 2 ./pagerank-parallel A.graph A.pr.txt	
	mpirun -np 2 ./pagerank-parallel B.graph B.pr.txt	
	mpirun -np 2 ./pagerank-parallel live-journal.graph live-journal.pr.txt
	mpirun -np 4 ./pagerank-parallel A.graph A.pr.txt	
	mpirun -np 4 ./pagerank-parallel B.graph B.pr.txt	
	mpirun -np 4 ./pagerank-parallel live-journal.graph live-journal.pr.txt
	mpirun -np 8 ./pagerank-parallel A.graph A.pr.txt	
	mpirun -np 8 ./pagerank-parallel B.graph B.pr.txt	
	mpirun -np 8 ./pagerank-parallel live-journal.graph live-journal.pr.txt
	mpirun -np 16 ./pagerank-parallel A.graph A.pr.txt	
	mpirun -np 16 ./pagerank-parallel B.graph B.pr.txt	
	mpirun -np 16 ./pagerank-parallel live-journal.graph live-journal.pr.txt
	mpirun -np 32 ./pagerank-parallel A.graph A.pr.txt	
	mpirun -np 32 ./pagerank-parallel B.graph B.pr.txt	
	mpirun -np 32 ./pagerank-parallel live-journal.graph live-journal.pr.txt
 
serial: pagerank.c pr_graph.c
	$(CC) -std=c99 -march=native -O3 -Wall -Wno-unused-result *.c -o pagerank -lm

clean:
	rm -f pagerank
