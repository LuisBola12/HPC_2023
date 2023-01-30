#include <mpi.h>
#include <stdio.h>

int main(int argc,char ** argv){
    int my_rank,num_of_process;
    float array_a[10000], array_b[10000];
    int error = -1;
    error = MPI_Init(&argc,&argv);
    MPI_Status status;  
    MPI_Comm_size(MPI_COMM_WORLD,&num_of_process);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    if (my_rank == 0){
        for(int i = 0; i < 10000; i++){
            array_a[i] = my_rank;
        }
        MPI_Send(array_a,10000,MPI_FLOAT,1,10,MPI_COMM_WORLD);
        MPI_Recv(array_b,10000,MPI_FLOAT,1,10,MPI_COMM_WORLD,&status);
        printf("I am task %d and I have received b[0] = %f",my_rank,array_b[0]);
    }else if(my_rank==1){
        for(int i = 0; i < 10000; i++){
            array_b[i] = my_rank;
        }
        MPI_Recv(array_a,10000,MPI_FLOAT,0,10,MPI_COMM_WORLD,&status);
        MPI_Send(array_b,10000,MPI_FLOAT,0,10,MPI_COMM_WORLD);
        printf("I am task %d and I have received a[0] = %f",my_rank,array_a[0]);
    }
    
    MPI_Finalize();

}