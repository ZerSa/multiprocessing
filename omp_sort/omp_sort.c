#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <math.h>

typedef struct scalar_ctx_t{
	int n;
	int m;
	int P;
	int N;
	double p;
} scalar_ctx_t;

void merge(int* input, int* output, int l1, int r1, int l2, int r2, int s){
    int i = 0, j = 0;
    while((l1 + i < r1 + 1) && (l2 + j < r2 + 1)){
        if (input[l1 + i] < input[l2 + j]){
            output[s + i + j] = input[l1 + i];
            i++;
        }
        else{
            output[s + i + j] = input[l2 + j];
            j++;
        }
    }
    
    while(l1 + i < r1 + 1){
        output[s + i + j] = input[l1 + i];
        i++;
    }
    
    while(l2 + j < r2 + 1){
        output[s + i + j] = input[l2 + j];
        j++;
    }
}

int binary_search (int *data, int x, int l_, int r_) {
    int l = l_;
    int r = r_;
    int h;
    if (l > r + 1)
        h = l;
    else
        h = r + 1;
        
    while (l < h) {
        int m = (l + h) / 2;
        if (x <= data[m]) {
            h = m;
        }
        else {
            l = m + 1;
        }
    }
    return h;
}
 

void swap (int* a, int* b) {
	int c = *a;
	*a = *b;
	*b = c;
}

void parallel_merge(int *input, int* output, int l1, int r1, int l2, int r2, int s){        
    if (r1 - l1 < r2 - l2){
        swap(&l1, &l2);
        swap(&r1, &r2);
    }
    
    if (r1 - l1 + 1 == 0)
        return;
    
    int m1 = (l1 + r1)/2;
    int m2 = binary_search(input, input[m1], l2, r2);
    int m3 = s + m1 - l1 + m2 - l2;
    
    output[m3] = input[m1];
    
    #pragma omp parallel
    {
        #pragma omp sections nowait
        {
            #pragma omp section
            {
                merge(input, output, l1, m1 - 1, l2, m2 - 1, s);
            }
            #pragma omp section
            {
                merge(input, output, m1 + 1, r1, m2, r2, m3 + 1);
            }
        }
    }
}




int comporator (const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
}

void omp_sort(int *input, int *output, int l, int r, int chunk_size, int s){
    int n = r - l + 1;
    if (r - l <= chunk_size){
        qsort(&input[l], n, sizeof(int), comporator);
        memcpy(&output[s], &input[l], n*sizeof(int));
        return;
    }
    
    int *tmp = (int *)malloc(n * sizeof(int)); 
    assert(tmp);
    int m = (l + r)/2;
     
    #pragma omp parallel
    {
        #pragma omp sections nowait
        {
            #pragma omp section
            {
                omp_sort(input, tmp, l, m, chunk_size, 0);
            }
            #pragma omp section
            {
                omp_sort(input, tmp, m + 1, r, chunk_size, m - l + 1);
            }
        }
    }
      
    parallel_merge(tmp, output, 0, m - l, m - l + 1, n - 1, s);
    free(tmp); 
}
    
int main (int argc, char **argv){
    if (argc != 4)
        return 0;
    
    int n = atoi(argv[1]);
	int m = atoi(argv[2]);
	int P = atoi(argv[3]);		

	int *data = calloc(n, sizeof(int));
	assert(data);
	int *sorted_data = calloc(n, sizeof(int)); 
	assert(sorted_data);
	
//	FILE *f_data = fopen("time_P.txt", "w");
	
	srand(time(NULL));
	for(int i = 0; i < n; ++i){
	    data[i] = rand()%1000;
	}
	
	omp_set_num_threads(P);
//	m = (int)(n*0.03);
	double start = omp_get_wtime();
	omp_sort(data, sorted_data, 0, n - 1, m, 0);  
    double end = omp_get_wtime();
    double omp_sort_time = end - start;
    
    FILE *f_stats = fopen("stats.txt", "w");
	FILE *f_data = fopen("data.txt", "w");
	
	for (int i = 0; i < n; ++i) {
		fprintf(f_data, "%d ", data[i]);
	}
	
	fprintf(f_data, "\n");
	
	for (int i = 0; i < n; ++i) {
		fprintf(f_data, "%d ", sorted_data[i]);
	}
	fprintf(f_data, "\n");
	
	fprintf(f_stats, "%fs %d %d %d \n", omp_sort_time, n , m, P);
	
	fclose(f_stats);
	//fprintf(f_data, "%f\n", omp_sort_time);
	
	
	/*
	double start = omp_get_wtime();
    qsort(data, n, sizeof(int), comporator);
	double end = omp_get_wtime();
	double quicksort_time = end - start;
	fprintf(f_data, "%f\n", quicksort_time);
	*/
	fclose(f_data);
	
	free(data);
	free(sorted_data);
	
	return 0;
}   
