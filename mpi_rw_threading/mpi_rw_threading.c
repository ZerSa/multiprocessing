#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <sys/time.h>
#include <mpi.h>

#define UP 1
#define DOWN 2
#define LEFT 3
#define RIGHT 4

#define MASTER 0
#define GROWTH 2
#define EXCHANGE_PERIOD 100


typedef struct ctx_type {
    int l;
    int a;
    int b;
    int n;
    int N;
    double p_u;
    double p_d;
    double p_l;
    double p_r;
    int rank;
    int size;
} ctx_t;

typedef struct point_t {
    int x;
    int y;
    int lifetime;
    int rank;
    
} point_t;

void push(point_t **data, int *size, int *capacity, point_t* elem){
    if (*size < *capacity){
        (*data)[*size] = *elem;
        (*size)++;
    }
    else{
        *data = realloc(*data, (*capacity) * GROWTH *sizeof(point_t));
        *capacity  = (*capacity) * GROWTH;
        (*data)[*size] = *elem;
        (*size)++;
    }
    return;
}

void pop(point_t **data, int *size, int index){
    if (*size != 0){ 
        (*data)[index] = (*data)[(*size) - 1];
        (*size)--;
    }
    return;
}


int create_seeds(int rank, int size){
    int *seeds;
    int seed;
    if (rank == MASTER){
        srand(time(NULL));
        seeds = malloc(size * sizeof(int));
        assert(seeds);
        for(int i = 0; i < size; ++i){
            seeds[i] = rand();
        }
    }
    
    MPI_Scatter(seeds, 1, MPI_INT, &seed, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    free(seeds);
    return seed;
}

int get_direction(void *ctx_){
    ctx_t *ctx = ctx_;
    double up = ctx->p_u * rand();
    double down = ctx->p_d * rand();
    double left = ctx->p_l * rand();
    double right = ctx->p_r * rand();
    if ((up >= down) && (up >= left) && (up >= right))
        return UP;
    if ((down >= up) && (down >= left) && (down >= right))
        return DOWN;
    if ((right >= up) && (right >= down) && (right >= left))
        return RIGHT;
    return LEFT;
} 

int get_next_process_rank(void *ctx_, int direction){
    ctx_t *ctx = ctx_;
    int y_rank = ctx->rank/ctx->a;
    int x_rank = ctx->rank % ctx->a;
    
    if(direction == UP){
        y_rank++;
        if (y_rank >= ctx->b)
            y_rank = 0;
    }
    if(direction == DOWN){
        y_rank--;
        if (y_rank < 0)
            y_rank = ctx->b - 1;
    }
    if(direction == RIGHT){
        x_rank++;
        if(x_rank >= ctx->a)
            x_rank = 0;
    }
    if(direction == LEFT){
        x_rank--;
        if(x_rank < 0)
            x_rank = ctx->a - 1;
    }
    
    return y_rank*ctx->a + x_rank;
}

void write_stats(void *ctx_, int *final_location, double start){
    ctx_t *ctx = ctx_;
    double end = MPI_Wtime();
    double delta = end - start;
    
    FILE *stats = fopen("stats.txt", "w");
    
    fprintf(stats, "%d %d %d %d %d %.2f %.2f %.2f %.2f %.4fs\n", ctx->l, ctx->a, ctx->b, 
            ctx->n, ctx->N, ctx->p_l, ctx->p_r, ctx->p_u, ctx->p_d, delta);
    
    fclose(stats);
    
}
    
    
void write_data(void *ctx_, point_t *points_completed, int size_completed){
    ctx_t *ctx = ctx_;
    MPI_File f;
    MPI_File_delete("data.bin", MPI_INFO_NULL);
    MPI_File_open(MPI_COMM_WORLD, "data.bin", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &f);
    
    int y_size = ctx->l;
    int x_size = ctx->l * ctx->size;
    
    int **result = calloc(y_size, sizeof(int *));
    assert(result);
    for(int i = 0; i < y_size; ++i){
        result[i] = calloc(x_size, sizeof(int));
        assert(result[i]);
    }
    
    for(int i = 0; i < size_completed; ++i){
        int y = points_completed[i].y;
        int x = points_completed[i].x;
        int initial_rank = points_completed[i].rank;
        result[y][x * ctx->size + initial_rank]++;
    }
    
    int a_global = ctx->rank % ctx->a;
    int b_global = ctx->rank / ctx->a;
    
    for(int i = 0; i < y_size; ++i){
        for(int j = 0; j < ctx->l; ++j){
            int number_bytes_line = ctx->l * ctx->size * ctx->a * sizeof(int);
            int number_bytes_square = number_bytes_line * b_global * ctx->l;
            int number_bytes_before_i = number_bytes_square + number_bytes_line * i;
            int number_bytes_first_line_elem = number_bytes_before_i + a_global*ctx->l *ctx->size*sizeof(int);
            int number_bytes_before_j = number_bytes_first_line_elem + j*ctx->size*sizeof(int);
            
            MPI_File_set_view(f, number_bytes_before_j, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
            MPI_File_write(f, &result[i][j*ctx->size], ctx->size, MPI_INT, MPI_STATUS_IGNORE);
        }
    }
    
    MPI_File_close(&f);
    
    for(int i = 0; i < y_size; ++i)
        free(result[i]);
        
    free(result);
}

    
    
    

void mpi_random_walk(void * ctx_){
    ctx_t *ctx = ctx_;
    
    int seed = create_seeds(ctx->rank, ctx->size);
    srand(seed);
    
    int local_points_number = ctx->N;
    int local_points_capacity = ctx->N;
    point_t *local_points = malloc(local_points_capacity * sizeof(point_t));
    assert(local_points);
    
    //создаем N точек в узле
    for(int i = 0; i < ctx->N; ++i){
        local_points[i].x = rand()%ctx->l;
        local_points[i].y = rand()%ctx->l;
        local_points[i].lifetime = 0;
        local_points[i].rank = ctx->rank;
    }
    
    // в эти массивы будем запихивать точки, которые перейдут в соседние узлы
    int size_send_left = 0;
    int capacity_send_left = ctx->N;
    point_t *points_send_left = malloc(capacity_send_left * sizeof(point_t));
    assert(points_send_left);
    
    
    int size_send_right = 0;
    int capacity_send_right = ctx->N;
    point_t *points_send_right = malloc(capacity_send_right * sizeof(point_t));
    assert(points_send_right);
    
    int size_send_down = 0;
    int capacity_send_down = ctx->N;
    point_t *points_send_down = malloc(capacity_send_down * sizeof(point_t));
    assert(points_send_down);
    
    int size_send_up = 0;
    int capacity_send_up = ctx->N;
    point_t *points_send_up = malloc(capacity_send_up * sizeof(point_t));
    assert(points_send_left);
    
    //а здесь будут точки, которые прошли свои N шагов и остались в этом узле
    int size_completed = 0;
    int capacity_completed = ctx->N;
    point_t *points_completed = malloc(capacity_completed * sizeof(point_t));
    assert(points_completed);
    
    int size_recv_up = 0;
    int size_recv_down = 0;
    int size_recv_left = 0;
    int size_recv_right = 0;
    
    
    int up_process_rank = get_next_process_rank(ctx, UP);
    int down_process_rank = get_next_process_rank(ctx, DOWN);
    int left_process_rank = get_next_process_rank(ctx, LEFT);
    int right_process_rank = get_next_process_rank(ctx, RIGHT);
    
    double start = MPI_Wtime();
    while(1){
        int current_ind = 0;
        while(current_ind < local_points_number){
            point_t *current_point = local_points + current_ind;
            int is_leaving = 0;
            
            for(int i = 0; i < EXCHANGE_PERIOD; i++){
                if (current_point->lifetime == ctx->n){
                    push(&points_completed, &size_completed, &capacity_completed, current_point);
                    pop(&local_points, &local_points_number, current_ind);
                    is_leaving = 1;
                    break;
                }
                
                current_point->lifetime++;
                int direction = get_direction(ctx);
                if(direction == UP)
                    current_point->y++;
                if (direction == DOWN)
                    current_point->y--;
                if(direction == LEFT)
                    current_point->x--;
                if(direction == RIGHT)
                    current_point->x++;
                    
                //проверяем, покидает ли точка пределы узла
                if (current_point->x >= ctx->l){
                    current_point->x = 0;
                    push(&points_send_right, &size_send_right, &capacity_send_right, current_point);
                    pop(&local_points, &local_points_number, current_ind);
                    is_leaving = 1;
                    break;
                }
                
                if (current_point->x < 0){
                    current_point->x = ctx->l - 1;
                    push(&points_send_left, &size_send_left, &capacity_send_left, current_point);
                    pop(&local_points, &local_points_number, current_ind);
                    is_leaving = 1;
                    break;
                }
                
                if (current_point->y >= ctx->l){
                    current_point->y = 0;
                    push(&points_send_up, &size_send_up, &capacity_send_up, current_point);
                    pop(&local_points, &local_points_number, current_ind);
                    is_leaving = 1;
                    break;
                }
                
                if (current_point->y < 0){
                    current_point->y = ctx->l - 1;
                    push(&points_send_down, &size_send_down, &capacity_send_down, current_point);
                    pop(&local_points, &local_points_number, current_ind);
                    is_leaving = 1;
                    break;
                }
            }
            if (is_leaving == 0)
                current_ind++;
        }
        
        //обмениваемся количество перешедших точек с соседними узлами
        MPI_Request* exchange_sizes = malloc(8 * sizeof(MPI_Request));
        assert(exchange_sizes);
        
        MPI_Isend(&size_send_left, 1, MPI_INT, left_process_rank, 0, MPI_COMM_WORLD, exchange_sizes + 0);
        MPI_Isend(&size_send_right, 1, MPI_INT, right_process_rank, 1, MPI_COMM_WORLD,    exchange_sizes + 1);
        MPI_Isend(&size_send_up, 1, MPI_INT, up_process_rank, 2, MPI_COMM_WORLD, exchange_sizes + 2);
        MPI_Isend(&size_send_down, 1, MPI_INT, down_process_rank, 3, MPI_COMM_WORLD, exchange_sizes + 3);
        
        MPI_Irecv(&size_recv_left, 1, MPI_INT, left_process_rank, 1, MPI_COMM_WORLD, exchange_sizes + 4);
        MPI_Irecv(&size_recv_right, 1, MPI_INT, right_process_rank, 0, MPI_COMM_WORLD,
exchange_sizes + 5);
        MPI_Irecv(&size_recv_up, 1, MPI_INT, up_process_rank, 3, MPI_COMM_WORLD, exchange_sizes + 6);
        MPI_Irecv(&size_recv_down, 1, MPI_INT, down_process_rank, 2, MPI_COMM_WORLD, exchange_sizes + 7);
                
                
        MPI_Waitall(8, exchange_sizes, MPI_STATUS_IGNORE);
        free(exchange_sizes);
    
        //обмениваемся самими перешедшими точками
        
        point_t* points_recv_left = malloc(size_recv_left* sizeof(point_t));
        assert(points_recv_left);

        point_t* points_recv_right = malloc(size_recv_right * sizeof(point_t));
        assert(points_recv_right);

        point_t* points_recv_up = malloc(size_recv_up * sizeof(point_t));
        assert(points_recv_up);

        point_t* points_recv_down = malloc(size_recv_down * sizeof(point_t));
        assert(points_recv_down);

        
        MPI_Request* exchange_points = malloc(8 * sizeof(MPI_Request));
        assert(exchange_points);
        
        MPI_Isend(points_send_left, size_send_left * sizeof(point_t), MPI_BYTE, left_process_rank, 0, MPI_COMM_WORLD, exchange_points + 0);
        MPI_Isend(points_send_right, size_send_right * sizeof(point_t), MPI_BYTE, right_process_rank, 1, MPI_COMM_WORLD, exchange_points + 1);
        MPI_Isend(points_send_up, size_send_up * sizeof(point_t), MPI_BYTE, up_process_rank, 2, MPI_COMM_WORLD, exchange_points + 2);
        MPI_Isend(points_send_down, size_send_down * sizeof(point_t), MPI_BYTE, down_process_rank, 3, MPI_COMM_WORLD, exchange_points + 3);
        
        MPI_Irecv(points_recv_left, size_recv_left * sizeof(point_t), MPI_BYTE,
left_process_rank, 1, MPI_COMM_WORLD, exchange_points + 4);
        MPI_Irecv(points_recv_right, size_recv_right * sizeof(point_t), MPI_BYTE,
right_process_rank, 0, MPI_COMM_WORLD, exchange_points + 5);
        MPI_Irecv(points_recv_up, size_recv_up * sizeof(point_t), MPI_BYTE,
up_process_rank, 3, MPI_COMM_WORLD, exchange_points + 6);
        MPI_Irecv(points_recv_down, size_recv_down * sizeof(point_t), MPI_BYTE, down_process_rank, 2, MPI_COMM_WORLD, exchange_points + 7);
        
        MPI_Waitall(8, exchange_points, MPI_STATUS_IGNORE);
        
        size_send_up = 0;
        size_send_down = 0;
        size_send_right = 0;
        size_send_left = 0;
        
        //теперь запихиваем полученные от других узлов точки в массив с текущими точками
        for (int i = 0; i < size_recv_up; ++i) {
            push(&local_points, &local_points_number,
                 &local_points_capacity, points_recv_up + i);
        }
        
        for (int i = 0; i < size_recv_down; ++i) {
            push(&local_points, &local_points_number,
                 &local_points_capacity, points_recv_down + i);
        }
        
        for (int i = 0; i < size_recv_left; ++i) {
            push(&local_points, &local_points_number,
                 &local_points_capacity, points_recv_left + i);
        }
        
        for (int i = 0; i < size_recv_right; ++i) {
            push(&local_points, &local_points_number,
                 &local_points_capacity, points_recv_right + i);
        }
        
        free(points_recv_up);
        free(points_recv_down);
        free(points_recv_left);
        free(points_recv_right);
        
        int is_finished = 0;
        int num_completed_points = 0;
        
        MPI_Reduce(&size_completed, &num_completed_points, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        
        if ((ctx->rank == 0) && (num_completed_points == ctx->size*ctx->N))
            is_finished = 1;
        
        MPI_Bcast(&is_finished, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
        
        if(is_finished == 1){
            int *final_location = calloc(ctx->size, sizeof(int));
            assert(final_location);
            
            MPI_Gather(&size_completed, 1, MPI_INT, final_location, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
            
            if(ctx->rank == 0){
                write_stats(ctx, final_location, start);
            }
            
            write_data(ctx, points_completed, size_completed);
            
            free(final_location);
            
            break;
            
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    free(points_send_up);
    free(points_send_down);
    free(points_send_left);
    free(points_send_right);
    
    free(local_points);
    free(points_completed);
}
        
    



int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    ctx_t ctx = {
        .l = atoi(argv[1]),
        .a = atoi(argv[2]),
        .b = atoi(argv[3]),
        .n = atoi(argv[4]),
        .N = atoi(argv[5]),
        .p_l = atof(argv[6]),
        .p_r = atof(argv[7]),
        .p_u = atof(argv[8]),
        .p_d = atof(argv[9]),
        .rank = rank,
        .size = size
    };
    
    mpi_random_walk(&ctx);
    MPI_Finalize();
    return 0;
}
