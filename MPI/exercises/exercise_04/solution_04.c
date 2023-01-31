#include <mpi.h>
#include <stdio.h>

int main(int argc,char ** argv){
    int my_rank,num_of_process;
    int error = -1;
    error = MPI_Init(&argc,&argv);
    MPI_Status status;  
    MPI_Comm_size(MPI_COMM_WORLD,&num_of_process);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    int my_left, my_right,my_sum, my_adding;
    my_sum = 0;
    my_left = (my_rank-1 < 0) ? num_of_process-1 : my_rank-1;
    my_right = (my_rank+1) % num_of_process;
    my_adding = my_rank;
    for(int i = 0; i < num_of_process; i++){
        MPI_Sendrecv(&my_adding,1,MPI_INT,my_left,10,&my_adding,1,MPI_INT,my_right,10,MPI_COMM_WORLD,&status);
        my_sum += my_adding;
    }
    printf("I am proc %d and sum = %d \n",my_rank,my_sum);
    MPI_Finalize();

}