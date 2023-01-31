#include <mpi.h>
#include <stdio.h>
#define FIXED_N 100

int main(int argc,char ** argv){
    int my_rank,num_of_process;
    int error = -1;
    
    error = MPI_Init(&argc,&argv);
    MPI_Status status;  
    
    MPI_Comm_size(MPI_COMM_WORLD,&num_of_process);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    
    int local_matrix_size = FIXED_N / num_of_process;
    int local_matrix[local_matrix_size + 2][FIXED_N];
    
    for(int i = 0;i<local_matrix_size;i++){
        for(int j =0;j<FIXED_N;j++){
            local_matrix[i][j]= my_rank;
        }
    }    
    if (my_rank > 0) {
        MPI_Send(local_matrix[1], FIXED_N, MPI_INT, my_rank-1, 10, MPI_COMM_WORLD);
        MPI_Recv(local_matrix[local_matrix_size+1], FIXED_N, MPI_INT, my_rank-1, 10, MPI_COMM_WORLD, &status);
    }
    if (my_rank < num_of_process-1) {
        MPI_Send(local_matrix[local_matrix_size], FIXED_N, MPI_INT, my_rank+1, 10, MPI_COMM_WORLD);
        MPI_Recv(local_matrix[0], FIXED_N, MPI_INT, my_rank+1, 10, MPI_COMM_WORLD, &status);
    }
    for(int i = 0;i<local_matrix_size;i++){
        for(int j =0;j<FIXED_N;j++){
            printf("%d ", local_matrix[i][j]);
    }
        printf("\n");
    }    
    MPI_Finalize();
}