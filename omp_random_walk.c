#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include <assert.h>
#include <time.h>

typedef struct scalar_ctx_t{
	int a;
	int b;
	int x;
	int N;
	double p;
} scalar_ctx_t;

void omp_random_walk(void *context, FILE *f){
	scalar_ctx_t *ctx = context;
	int walking_time = 0;
	int b_end_counter = 0; 
	
	int *seed = (int *)malloc(ctx->N * sizeof(int));
	assert(seed);
	int s = (int)(time(NULL));
	for(int i = 0; i < ctx->N; ++i){
		seed[i] = rand_r(&s);
	}
	
	struct timeval start, end;
	assert(gettimeofday(&start, NULL) == 0);
	
	#pragma omp parallel for reduction(+: walking_time)
	for (int i = 0; i < ctx->N; ++i){
		int local_x = ctx->x;
		int local_seed = seed[i];
		while( local_x != ctx->a && local_x != ctx->b){
			double p = (double)rand_r(&local_seed)/ (RAND_MAX);
			walking_time++;
			if (p <= ctx->p){
				local_x++;
			}
			else{
				local_x--;
			}
		}
		
		if (local_x == ctx->b){
			#pragma omp atomic
			b_end_counter++;
		}
	}
	
	assert(gettimeofday(&end, NULL) == 0);
	double delta = ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec -
		start.tv_usec) / 1.e6;
	double probability_b = (double)b_end_counter/ctx->N;
	double mean_walking_time = (double)walking_time/ctx->N;
	fprintf(f, "%f %f %fs %d %d %d %d %f \n", probability_b, mean_walking_time, delta, ctx->a, ctx->b, ctx->x, ctx->N, ctx->p);
	free(seed);
}

int main(int argc, char **argv){
	if (argc != 7){
		return 0;
	}
	
	int a = atoi(argv[1]);
	int b = atoi(argv[2]);
	int x = atoi(argv[3]);
	int N = atoi(argv[4]);
	double p = atof(argv[5]);
	int P = atoi(argv[6]);

	scalar_ctx_t ctx = {
		.a = a,
		.b = b,
		.x = x,
		.N = N,
		.p = p,
	};
	
	omp_set_num_threads(P);
	FILE *f = fopen("stats.txt", "w");
	if (f != NULL){
		omp_random_walk(&ctx, f);
	}
	
	fclose(f);
	
	
	return 0;
}
	 
