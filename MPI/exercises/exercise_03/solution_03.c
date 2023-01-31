#include <mpi.h>
#include <stdio.h>

int main(int argc,char ** argv){
    int my_rank,num_of_process;
    float my_array[10000], partner_array[10000];
    int error = -1;
    error = MPI_Init(&argc,&argv);
    MPI_Status status;  
    MPI_Comm_size(MPI_COMM_WORLD,&num_of_process);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    int my_left, my_right;
    my_left = (my_rank-1 < 0) ? num_of_process-1 : my_rank-1;
    my_right = (my_rank+1) % num_of_process;
    for(int i = 0; i < 10000; i++){
        my_array[i] = my_rank;
    }
    MPI_Sendrecv(&my_array,10000,MPI_FLOAT,my_right,10,&partner_array,10000,MPI_FLOAT,my_left,10,MPI_COMM_WORLD,&status);
    printf("I am task %d and I have received partner[0] = %f \n",my_rank,partner_array[0]);
    MPI_Finalize();

}