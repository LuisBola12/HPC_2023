#include <mpi.h>
#include <stdio.h>

int main(int argc,char ** argv){
    int my_rank,num_of_process;
    int error = -1;
    error = MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&num_of_process);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    printf("Hello world, I am proc %d of total %d \n",my_rank,num_of_process);
}