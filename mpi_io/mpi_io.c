#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <assert.h>

typedef struct scalar_ctx_t {
	int l;
	int a;
	int b;
	int N;
	int* d;
} scalar_ctx_t;


void mpi_io(void *context){
    scalar_ctx_t *ctx = context;
    int seed = time(NULL);
    int comm_rank;
    struct timespec start, end;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
    if (!comm_rank){
        int  comm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	
	if (comm_size != ctx->a * ctx->b) {
	    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}
        
        int *seeds = malloc(comm_size * sizeof(int));
        assert(seeds);
        
	int i;
        for(i = 0; i < comm_size; ++i){
            seeds[i] = seed + i;
        }
        
        MPI_Scatter(seeds, 1, MPI_INT, &seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        free(seeds);
        
        clock_gettime(CLOCK_MONOTONIC, &start);
    }
    else{
        MPI_Scatter(NULL, 1, MPI_INT, &seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    
    srand(seed);
    int i;
    for(i = 0; i < ctx->N; ++i){
        int x = rand() % ctx->l;
        int y = rand() % ctx->l;
        int r = rand() % ctx->a*ctx->b;
        ctx->d[(y * ctx->l + x) *ctx->a * ctx->b + r] += 1;
    }
    
    MPI_File data;
    MPI_Datatype contiguous, view;
    
    MPI_File_open(MPI_COMM_WORLD, "data.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &data);
    MPI_File_set_size(data, 0);
    
    MPI_Type_contiguous(ctx->l * ctx->a *ctx->b, MPI_INT, &contiguous);
    MPI_Type_create_resized(contiguous, 0, ctx->a*ctx->l*ctx->a*ctx->b*sizeof(int), &view);
    MPI_Type_commit(&view);
    MPI_File_set_view(data, ((comm_rank / ctx->a) * ctx->a * ctx->l + comm_rank % ctx->a) * ctx->l * ctx->a*ctx->b * sizeof(int), MPI_INT, view, "native", MPI_INFO_NULL);
	MPI_File_write_all(data, ctx->d, ctx->l * ctx->l * ctx->a*ctx->b, MPI_INT, MPI_STATUS_IGNORE);

	MPI_File_close(&data);
	free(ctx->d);
	
	if(!comm_rank){
	    clock_gettime(CLOCK_MONOTONIC, &end);
	    double delta = end.tv_sec - start.tv_sec;
	    delta += (end.tv_nsec - start.tv_nsec) / 1000000000.0;
	    
	    FILE *stats = fopen("stats.txt", "w");
	    fprintf(stats, "%d %d %d %d %fs\n", ctx->l, ctx->a, ctx->b, ctx->N, delta);
	    fclose(stats);
	}
}
      

int main(int argc, char **argv){
    MPI_Init(&argc, &argv);
    
    if (argc != 5)
        return 0;
    
    scalar_ctx_t ctx = {
        .l = atoi(argv[1]),
		.a = atoi(argv[2]),
		.b = atoi(argv[3]),
		.N = atoi(argv[4]),
		.d = 0
	}; 
	
	ctx.d = calloc(ctx.l * ctx.l * ctx.a * ctx.b, sizeof(int));
	mpi_io(&ctx);
	
	MPI_Finalize();
	return 0;
}
