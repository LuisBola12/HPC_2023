#include <mpi.h>
#include <stdio.h>
#define FIXED_N 10

int main(int argc,char ** argv){
    int my_rank,num_of_process;
    int error = -1;
    
    error = MPI_Init(&argc,&argv);
    MPI_Status status;  
    
    //MPI_Comm_size(MPI_COMM_WORLD,&num_of_process);
    num_of_process = 4;
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    
    int matrix[FIXED_N][FIXED_N];
    
    for(int i = 0;i<FIXED_N;i++){
        for(int j =0;j<FIXED_N;j++){
            matrix[i][j] = -1;
        }
    }
    int reminder = FIXED_N%num_of_process;
    int work_to_do = FIXED_N/num_of_process;
    int work[FIXED_N];
    if(my_rank == 0){
        int current_index = 0;
        if(reminder > 0){
            for(int i = 0; i < num_of_process;i++){
                if(reminder > 0){
                    for(int j = 0; j < work_to_do+1;j++){
                        work[current_index]=i;
                        current_index++;
                    }
                    reminder--;
                    
                }else{
                   for(int j = 0; j < work_to_do;j++){
                        work[current_index]=i;
                        current_index++;
                    } 
                }
            }
        }else{
            for(int i = 0; i < num_of_process;i++){            
               for(int j = 0; j < work_to_do;j++){
                    work[current_index]=i;
                    current_index++;
                } 
            }
        }     
    }
    MPI_Bcast(work,FIXED_N,MPI_INT,10,MPI_COMM_WORLD);
    int n_work = 0
    for(int i = 0; i < FIXED_N;i++){
        if(work[i]==my_rank){
            n_work++;
        }
    }
    int my_work[n_work];
    int current_index = 0;
    for(int i = 0; i < FIXED_N;i++){
        if(work[i]==my_rank){
            my_work[current_index] = i;
            current_index++;
        }
    }
    for(int i = 0; i < n_work; i++){
        matrix[diagonal][diagonal] = my_work[i]
    }
    if(my_rank ==0){

    }else{
        
    }
    MPI_Finalize();
}