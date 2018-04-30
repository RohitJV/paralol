
pagerank: pagerank-parallel.c pr_graph.c
	mpicc pagerank-parallel.c -o pagerank-parallel	
	mpirun -np 16 ./pagerank-parallel B-small.graph A-small-parallel.txt

clean:
	rm -f pagerank
