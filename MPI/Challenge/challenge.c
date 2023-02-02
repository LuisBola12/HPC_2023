/*
 *    It evolves the equation:
 *                            u,t + u,x + u,y = 0
 *    Using a Lax scheme.
 *    The initial data is a cruddy gaussian.
 *    Boundaries are flat: copying the value of the neighbour
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define NX 100
#define NY 100
#define LX 2.0
#define LY 2.0

#define sigma 0.1
#define tmax 100


/* 
 * conversions from discrete to real coordinates
 */
float ix2x(int ix){
    return ((ix-1) - (NX-1) / 2.0)*LX/(NX-1);
}

float iy2y(int iy){ 
    return ((iy-1) - (NY-1) / 2.0)*LY/(NY-1);
}


/* Function for evaluating the results that calculates the sum of the array values */
float summa(int nx, int ny, float* val){
float summma=0.0;
int ix,iy;
    for(iy=1;iy<=ny;++iy)
	for(ix=1;ix<=nx;++ix){
          summma+=val[((nx+2)*iy)+ix];
 } 
  return(summma);
}


/* 
 * initialize the system with a gaussian temperature distribution 
 */
int init_transport(float *temp){ 
    int ix,iy;
    float x,y;

    for(iy=1;iy<=NY;++iy)
	for(ix=1;ix<=NX;++ix){
	    x=ix2x(ix);
	    y=iy2y(iy);
	    temp[((NX+2)*iy)+ix] = tmax*exp((-((x*x)+(y*y)))/(2.0*(sigma*sigma)));
	}
    return 0;
}

/*
 * save the temperature distribution
 * the ascii format is suitable for splot gnuplot function
 */
int save_gnuplot(char *filename, float *temp,int my_rank,int num_of_rows,int num_procs) {
    int ix,iy;
    float * row;
    FILE *fp;
    MPI_Status status;
    if(my_rank == 0){
        fp = fopen(filename, "w"); 
        for(iy=1;iy<=num_of_rows;++iy){
	        for(ix=1;ix<=NX;++ix){
	            fprintf(fp, "\t%f\t%f\t%g\n", ix2x(ix), iy2y(iy), temp[((NX+2)*iy)+ix]);
            }
            fprintf(fp, "\n");
        }
        row = (float *) malloc ((NX+2)*(num_of_rows+2)*sizeof(float));
        for(int i =1;i<num_procs;i++){
            MPI_Recv(row,((NX+2)*(NY+2)),MPI_REAL,i,1,MPI_COMM_WORLD,&status);
            for(iy=1;iy<=num_of_rows;++iy){
	            for(ix=1;ix<=NX;++ix){
	                fprintf(fp, "\t%f\t%f\t%g\n", ix2x(ix), iy2y(iy), row[((NX+2)*iy)+ix]);
            }
            fprintf(fp, "\n");
        }
        }
        fclose(fp);
    }else{
        MPI_Send(temp,((NX+2)*(NY+2)),MPI_REAL,0,1,MPI_COMM_WORLD);
    }  
    
    return 0;
}

int update_boundaries_FLAT(float *temp,int num_proces,int my_rank,int num_rows,int rank_up, int rank_bhd){
    int ix=0, iy=0;
    //QUEDA IGUAL POR QUE SON LAS COLUMNAS
    for(iy=1;iy<=num_rows;++iy){
	    temp[(NX+2)*iy] = temp[((NX+2)*iy)+1];
	    temp[((NX+2)*iy)+(NX+1)] = temp[((NX+2)*iy)+NX];
    }
    //SOLO PARA LOS EXTREMOS
    if(my_rank == 0){
        //LLENA SOLO LA PRIMERA
        for(ix=0;ix<=num_rows+1;++ix){
	        temp[ix] = temp[(NX+2)+ix];
        }
    }else if(my_rank == num_proces-1){
        //LA ULTIMA FILA
        for(ix=0;ix<=num_rows+1;++ix){
	        temp[((NX+2)*(num_rows+1))+ix] = temp[((NX+2)*num_rows)+ix];
        }
    }
    MPI_Status  status1, status2;
    MPI_Sendrecv(&temp[(NX+2)+1],NX,MPI_REAL,rank_up,10,
      &temp[((NX+2)*(num_rows+1))+1],NX,MPI_REAL,rank_bhd,10,
      MPI_COMM_WORLD,&status1);

    MPI_Sendrecv(&temp[((NX+2)*num_rows)+1],NX,MPI_REAL,rank_bhd,
      11,&temp[1],NX,MPI_REAL,rank_up,11,MPI_COMM_WORLD,&status2);
    return 0;
}


int evolve(float dtfact, float *temp, float *temp_new,int num_of_rows){
    float dx, dt;
    int ix, iy;
    float temp0;

    dx = 2*LX/NX;
    dt = dtfact*dx/sqrt(3.0);
    for(iy=1;iy<=num_of_rows;++iy)
	for(ix=1;ix<=NX;++ix){
	    temp0 = temp[((NX+2)*iy)+ix];
	    temp_new[((NX+2)*iy)+ix] = temp0-0.5*dt*(temp[((NX+2)*(iy+1))+ix]
            -temp[((NX+2)*(iy-1))+ix]+temp[((NX+2)*iy)+(ix+1)]-temp[((NX+2)*iy)+(ix-1)])/dx;
	}

    for(iy=0;iy<=num_of_rows+1;++iy)
	for(ix=0;ix<=NX+1;++ix)
	    temp[((NX+2)*iy)+ix] = temp_new[((NX+2)*iy)+ix];
    

    return 0;
}

int main(int argc, char* argv[]){
    int i=0, nRow=NX+2, nCol=NY+2;
    int num_process,my_rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&my_rank);
    MPI_Comm_rank(MPI_COMM_WORLD,&num_process);

    int rank_up,rank_bhd;
    rank_up = my_rank + 1;
    rank_bhd = my_rank - 1;


    if(my_rank == 0){
        rank_bhd = MPI_PROC_NULL;
    }
    if(my_rank == num_process-1){
        rank_up = MPI_PROC_NULL;
    }
    
    int number_of_rows = NY/num_process;
    if (my_rank < NY % num_process){
        number_of_rows++;
    }

    float *temp, *temp_new;
    temp = (float *) malloc (nRow*nCol*sizeof(float));
    temp_new = (float *) malloc (nRow*nCol*sizeof(float));

    init_transport(temp);
    update_boundaries_FLAT(temp,num_process,my_rank,number_of_rows,rank_up,rank_bhd);

    float before=summa(NX, number_of_rows, temp);
    printf(" sum temp before: %f\n",before);

    save_gnuplot("transport.dat", temp,my_rank,number_of_rows,num_process);

    for(i=1;i<=500;++i) {
        evolve(0.1, temp, temp_new,number_of_rows);
        update_boundaries_FLAT(temp,num_process,my_rank,number_of_rows,rank_up,rank_bhd);
    }

    float after=summa(NX, number_of_rows, temp);
    printf(" sum temp after: %f\n",after);
    save_gnuplot("transport_end.dat", temp,my_rank,number_of_rows,num_process);
 
  
    free(temp);
    free(temp_new);    
    return 0;
}