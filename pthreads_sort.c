#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>

int n;
int m;
int P;

typedef struct pthread_data_t {
    int l1; 
    int r1;
    int l2; 
    int r2;
    int* data_location;
} pthread_data_t;

typedef struct chunks_data_t {
    int start;
    int chunks_number;
    int* chunks_location;
} chunks_data_t; 
    
    
int comparator (const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
}


void simple_merge(int* input, int l1, int r1, int l2, int r2){
    int i = 0, j = 0;
    int n1 = r1 - l1 + 1;
    int n2 = r2 - l2 + 1;
    
    int *tmp1 = (int*)malloc(n1*sizeof(int));
    assert(tmp1);
    memcpy(tmp1, &input[l1], n1*sizeof(int));
    
    int *tmp2 = (int*)malloc(n2*sizeof(int));
    assert(tmp2);
    memcpy(tmp2, &input[l2], n2*sizeof(int));
    
    while((i < n1) && (j < n2)){
        if(tmp1[i] < tmp2[j]){
            input[l1 + i + j] = tmp1[i];
            i++;
        }
        else{
            input[l1 + i + j] = tmp2[j];
            j++;
        }
    }
    
    while(i < n1){
        input[l1 + i + j] = tmp1[i];
        i++;
    }
    
    while(j < n2){
        input[l1 + i + j] = tmp2[j];
        j++;
    }
    
    free(tmp1);
    free(tmp2);
}
    
//функция для потока для сортировки чанков
void* thread_chunks_sort(void *chunks_data_){
    chunks_data_t *chunks_data = (chunks_data_t *)chunks_data_;
    int start = chunks_data->start;
    int chunks_number = chunks_data->chunks_number;
    int* chunks_location = chunks_data->chunks_location;
    
    int current_ind = start;
    
    for(int i = 0; i < chunks_number; ++i){
        qsort(&chunks_location[current_ind], m, sizeof(int), comparator);
        current_ind+= m;
    }
    return NULL;
}

//функция для потока для слияния
void * thread_simple_merge(void* data_){
    pthread_data_t *data = (pthread_data_t *)data_;
    simple_merge(data->data_location, data->l1, data->r1, data->l2, data->r2);
    return NULL;
}

void parallel_sort_chunks(int * input, int n){
    pthread_t* threads = (pthread_t *)malloc(P * sizeof(pthread_t));
    assert(threads);
    chunks_data_t *thread_data = malloc(P * sizeof(chunks_data_t));
    assert(thread_data);
    
    int chunks_number_for_thread = n / m / P;
    int current_ind = 0;
    
    // создаем потоки, которые сортируют чанки
    for(int i = 0; i < P; ++i){
        thread_data[i].start = current_ind;
        thread_data[i].chunks_number = chunks_number_for_thread;
        thread_data[i].chunks_location = input;
            
        pthread_create(&(threads[i]), NULL, thread_chunks_sort, &thread_data[i]);
        current_ind+= m* chunks_number_for_thread;
    }
         
    for(int i = 0; i < P; ++i){
        pthread_join(threads[i], NULL);
    }
    
    qsort(&input[current_ind], n - current_ind, sizeof(int), comparator);
    
    free(threads);
    free(thread_data);
}



void pthread_sort(int *input, int n){
    parallel_sort_chunks(input, n);
    
    pthread_t *threads = malloc(P * sizeof(pthread_t));
    assert(threads);
    
    pthread_data_t *thread_data = malloc(P * sizeof(pthread_data_t));
    assert(thread_data);
    
    //слияние
    for(int current_size = m; current_size < n; current_size*=2){
        int current_ind = 0;
        int num_threads = P;
        int l1, l2, r1, r2;
        
        while(current_ind + current_size < n - 1){
            for(int i = 0; i < P; ++i){
                l1 = current_ind;
                r1 = current_ind + current_size -1;
                if (r1 >= n){
                    num_threads = i;
                    break;
                }
                
                l2 = r1 + 1;
                r2 = l2 + current_size - 1;
                if (r2 > n - 1)
                    r2 = n - 1;
                    
                current_ind = r2 + 1;
                
                thread_data[i].l1 = l1;
                thread_data[i].r1 = r1;
                thread_data[i].l2 = l2;
                thread_data[i].r2 = r2;
                thread_data[i].data_location = input;
                
                pthread_create(&(threads[i]), NULL, thread_simple_merge, &thread_data[i]);
            }
            
            for(int i = 0; i < num_threads; ++i){
                pthread_join(threads[i], NULL);
            }
        }
    }
    
    free(threads);
    free(thread_data);
}
                
        
    
    

int main (int argc, char **argv){
    if (argc != 4)
        return 0;
    
    n = atoi(argv[1]);
	m = atoi(argv[2]);
	P = atoi(argv[3]);		

	int *data = calloc(n, sizeof(int));
	assert(data);
	int *sorted_data = calloc(n, sizeof(int)); 
	assert(sorted_data);
	
//    FILE *f_stats = fopen("time_P.txt", "w");
	
	srand(time(NULL));
	for(int i = 0; i < n; ++i){
	    data[i] = rand()%1000;
	}
	
	FILE *f_stats = fopen("stat.txt", "w");
	FILE *f_data = fopen("data.txt", "w");
//	for (int p = 1; p <= 16; p = p*2){
//	P = p;
//	m = (int)(n*0.03);
	for (int i = 0; i < n; ++i) {
		fprintf(f_data, "%d ", data[i]);
	}
	fprintf(f_data, "\n");
	
	struct timeval start, end;
	assert(gettimeofday(&start, NULL) == 0);
	
	pthread_sort(data, n);  
	
	assert(gettimeofday(&end, NULL) == 0);
	double pthread_time = ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    
	
	for (int i = 0; i < n; ++i) {
		fprintf(f_data, "%d ", data[i]);
	}
	
	fprintf(f_data, "\n");
	
	fprintf(f_stats, "%fs %d %d %d \n", pthread_time, n , m, P);
//	fprintf(f_stats, "%f\n", pthread_time);
//	}
	
	
	
	fclose(f_data);
	fclose(f_stats);
	
	free(data);
	
	return 0;
}   
