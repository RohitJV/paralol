pagerank-parallel: pagerank-parallel.c pr_graph.c
	mpicc -lm pagerank-parallel.c -o pagerank-parallel
	mpirun -np 16 ./pagerank-parallel A.graph A.pr.txt

pagerank-parallel-improved: pagerank-parallel-improved.c pr_graph.c
	mpicc -std=c99 -march=native -O3 -lm pagerank-parallel-improved.c -o pagerank-parallel-improved
	mpirun -np 1 ./pagerank-parallel-improved A.graph A.pr.txt
	diff A.pr.txt A-serial.pr.txt
	mpirun -np 2 ./pagerank-parallel-improved A.graph A.pr.txt
	diff A.pr.txt A-serial.pr.txt
	mpirun -np 4 ./pagerank-parallel-improved A.graph A.pr.txt
	diff A.pr.txt A-serial.pr.txt
	mpirun -np 8 ./pagerank-parallel-improved A.graph A.pr.txt
	diff A.pr.txt A-serial.pr.txt
	mpirun -np 16 ./pagerank-parallel-improved A.graph A.pr.txt
	diff A.pr.txt A-serial.pr.txt
	mpirun -np 32 ./pagerank-parallel-improved A.graph A.pr.txt
	diff A.pr.txt A-serial.pr.txt
	mpicc -std=c99 -march=native -O3 -lm pagerank-parallel-improved.c -o pagerank-parallel-improved
	mpirun -np 1 ./pagerank-parallel-improved B.graph B.pr.txt
	diff B.pr.txt B-serial.pr.txt
	mpirun -np 2 ./pagerank-parallel-improved B.graph B.pr.txt
	diff B.pr.txt B-serial.pr.txt
	mpirun -np 4 ./pagerank-parallel-improved B.graph B.pr.txt
	diff B.pr.txt B-serial.pr.txt
	mpirun -np 8 ./pagerank-parallel-improved B.graph B.pr.txt
	diff B.pr.txt B-serial.pr.txt
	mpirun -np 16 ./pagerank-parallel-improved B.graph B.pr.txt
	diff B.pr.txt B-serial.pr.txt
	mpirun -np 32 ./pagerank-parallel-improved B.graph B.pr.txt
	diff B.pr.txt B-serial.pr.txt
	mpicc -std=c99 -march=native -O3 -lm pagerank-parallel-improved.c -o pagerank-parallel-improved
	mpirun -np 1 ./pagerank-parallel-improved live-journal.graph live-journal.pr.txt
	diff live-journal.pr.txt live-journal-serial.pr.txt
	mpirun -np 2 ./pagerank-parallel-improved live-journal.graph live-journal.pr.txt
	diff live-journal.pr.txt live-journal-serial.pr.txt
	mpirun -np 4 ./pagerank-parallel-improved live-journal.graph live-journal.pr.txt
	diff live-journal.pr.txt live-journal-serial.pr.txt
	mpirun -np 8 ./pagerank-parallel-improved live-journal.graph live-journal.pr.txt
	diff live-journal.pr.txt live-journal-serial.pr.txt
	mpirun -np 16 ./pagerank-parallel-improved live-journal.graph live-journal.pr.txt
	diff live-journal.pr.txt live-journal-serial.pr.txt
	mpirun -np 32 ./pagerank-parallel-improved live-journal.graph live-journal.pr.txt
	diff live-journal.pr.txt live-journal-serial.pr.txt

pagerank: pagerank.c pr_graph.c
	$(CC) -std=c99 -march=native -O3 pagerank.c pr_graph.c -o pagerank -lm
	./pagerank A.graph A-serial.pr.txt

clean:
	rm -f pagerank
	rm -f pagerank-parallel
	rm -f pagerank-parallel-improved
