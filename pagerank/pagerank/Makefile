
pagerank: pagerank.c pr_graph.c
	$(CC) -std=c99 -march=native -O3 -Wall -Wno-unused-result *.c -o pagerank -lm	
	./pagerank ../pagerank-parallel/A.graph A.pr.txt	
	./pagerank ../pagerank-parallel/B.graph B.pr.txt
	./pagerank ../pagerank-parallel/live-journal.graph live-journal.pr.txt

clean:
	rm -f pagerank
