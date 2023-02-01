//  
//   Inspired by Advanced MPI (Archer Training Material,EPCC)
//  
//   http://www.archer.ac.uk/training/course-material/2016/09/160929_AdvMPI_EPCC/index.php
//   
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

// Boundary values 
#define TOP 1.0
#define BOTTOM 10.0
#define LEFT 1.0
#define RIGHT 1.0
//
// The maximum number of iterations
#define MAX_ITERATIONS 5000000
// The convergence to terminate at
#define CONVERGENCE_ACCURACY 1e-4
// How often to report the norm
#define REPORT_NORM_PERIOD 1000


int nx, ny, ny2;

void initialise(double**, double**, int);
double* allocate_matrix_as_array(int nrows, int ncols);
double** allocate_matrix(int nrows, int ncols, double* arr_A);

int main(int argc, char * argv[]) {	

	if (argc != 3) {
		fprintf(stderr, "You must provide the size in X and size in Y as arguments to this code\n");
		return -1;
	}
	nx=atoi(argv[1]);
	ny=atoi(argv[2]);
        ny2=ny+2;

	printf("Solving to accuracy of %.0e, global system size is x=%d y=%d\n", CONVERGENCE_ACCURACY, nx, ny);
	int local_nx=nx; //This has to be changed in case of MPI

	double * grid1d =  allocate_matrix_as_array(local_nx+2,ny+2);
	double * grid_new1d = allocate_matrix_as_array(local_nx+2,ny+2);

        double **grid = allocate_matrix(local_nx+2,ny+2,grid1d);
        double **grid_new = allocate_matrix(local_nx+2,ny+2,grid_new1d);
	double start_time;

	initialise(grid, grid_new, local_nx);

	double rnorm=0.0, bnorm=0.0, norm, tmpnorm=0.0;

        // Initial norm factor
	int i,j,k;
	for (i=1;i<=local_nx;i++) {
		for (j=1;j<ny+1;j++) {		
			tmpnorm=tmpnorm+pow(grid[i][j]*4-grid[i][j-1]-grid[i][j+1]-grid[i-1][j]-grid[i+1][j],2);
		}
	}
	bnorm=sqrt(tmpnorm);

        //main loop
	for (k=0;k<MAX_ITERATIONS;k++) {
		tmpnorm=0.0;
    		//check convergence
		for (i=1;i<=local_nx;i++) {
			for (j=1;j<ny+1;j++) {
				tmpnorm=tmpnorm+pow(grid[i][j]*4-grid[i][j-1]-grid[i][j+1]-grid[i-1][j]-grid[i+1][j],2);
			}
		}
		norm=sqrt(tmpnorm)/bnorm;
		//Update grid
		if (norm < CONVERGENCE_ACCURACY) break;		
		for (i=1;i<=local_nx;i++) {
			for (j=1;j<ny+1;j++) {
				grid_new[i][j]=0.25*(grid[i][j-1]+grid[i][j+1]+grid[i-1][j]+grid[i+1][j]);
			}
		}
		memcpy(grid1d, grid_new1d, sizeof(double) * (local_nx + 2) * ny2);

		if (k % REPORT_NORM_PERIOD == 0) printf("Iteration= %d Relative Norm=%e\n", k, norm);
	}
	printf("\nTerminated on %d iterations, Relative Norm=%e\n", k, norm);
	free(grid1d);
	free(grid_new1d);
	free(grid);
	free(grid_new);
	return 0;
}

/**
 * Initialises the arrays, such that grid contains the boundary conditions at the start and end points and all other
 * points are zero. grid_new is set to equal grid
 */
void initialise(double ** grid, double ** grid_new, int local_nx) {
	int i, j;
	clock_t start, end;
	double elapsed;
	start = clock();
	#pragma omp parallel for 
    for (j=0;j<=local_nx+1;j++) {
        grid_new[j][0]=grid[j][0]=LEFT;
    }
	#pragma omp parallel for 
    for (j=0;j<=local_nx+1;j++) {
    	grid_new[j][ny2-1]=grid[j][ny2-1]=RIGHT;
    }
	#pragma omp parallel for 
	for (j=0;j<=ny+1;j++) {
		grid_new[0][j]=grid[0][j]=TOP;
	}
	#pragma omp parallel for 
	for (j=0;j<=ny+1;j++) {
		grid_new[local_nx+1][j]=grid[local_nx+1][j]=BOTTOM;		
	}
	#pragma omp parallel for collapse(2)
	for (i=1;i<=local_nx;i++) {
		for (j=1;j<ny+1;j++) {
			grid_new[i][j]=grid[i][j]=0;
		}
	}
	end = clock();
	elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
  	printf("Time elapsed: %f\n", elapsed);
}


/* Allocate a double matrix with one malloc */
double* allocate_matrix_as_array(int nrows, int ncols) {
  double *arr_A;
  /* Allocate enough raw space */
  arr_A = (double*) malloc(nrows*ncols*sizeof(double));
  return arr_A;
}

/* Easy way to access a 1d dynamic array */
double** allocate_matrix(int nrows, int ncols, double* arr_A) {
  double **A;
  int ii;
  /* Prepare pointers for each matrix row */
  A = (double**) malloc(sizeof(double)*nrows);
  /* Initialize the pointers */
  for (ii=0; ii<nrows; ++ii) {
    A[ii] = &(arr_A[ii*ncols]);
  }
  return A;
}

