pagerank-parallel: pagerank-parallel.c pr_graph.c
	mpicc -lm pagerank-parallel.c -o pagerank-parallel
	mpirun -np 16 ./pagerank-parallel A.graph A.pr.txt

pagerank-parallel-improved: pagerank-parallel-improved.c pr_graph.c
	mpicc -lm pagerank-parallel-improved.c -o pagerank-parallel-improved
	mpirun -np 16 ./pagerank-parallel-improved A.graph A.pr.txt		

pagerank: pagerank.c pr_graph.c
	$(CC) -std=c99 -march=native -O3 -Wall -Wno-unused-result pagerank.c pr_graph.c -o pagerank -lm
	./pagerank A.graph A-serial.pr.txt

clean:
	rm -f pagerank
	rm -f pagerank-parallel
	rm -f pagerank-parallel-improved
