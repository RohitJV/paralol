default:
	gcc -pthread lc_pthreads.c -o lc_pthreads
	gcc -fopenmp lc_openmp.c -o lc_openmp
	./lc_pthreads data_large.csv data_large.label 10 1
	./lc_pthreads data_large.csv data_large.label 10 2
	./lc_pthreads data_large.csv data_large.label 10 4
	./lc_pthreads data_large.csv data_large.label 10 8
	./lc_pthreads data_large.csv data_large.label 10 16
	./lc_openmp data_large.csv data_large.label 10 1
	./lc_openmp data_large.csv data_large.label 10 2
	./lc_openmp data_large.csv data_large.label 10 4
	./lc_openmp data_large.csv data_large.label 10 8
	./lc_openmp data_large.csv data_large.label 10 16
