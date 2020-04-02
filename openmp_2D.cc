#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include "mpi.h"
#define M 3360
#define W 29

using namespace std;

/* window matrix as global vabriable */
float h[W][W];

/** main function ***/
int main(int argc, char* argv[]){

   
 double start, end;             /* timestamps */
 int i, j, q;    
 MPI_Status status;

/* initialize the window matrix */
  for ( int i = 0; i < W; i++){
    for (int j = 0; j < W; j++){
         h[i][j] = 1.0/(W*W);
  }
 }



 MPI_Init(&argc, &argv);

 int my_rank, p;

 MPI_Comm_size(MPI_COMM_WORLD, &p);

 MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

 int N = (M/p) + (W - 1);       /* row size of comm buffer */
 int Nm = (M/p) + ((W - 1)/2);  /* size of partition for last process */


/* rank 0 process generates image matrix f */
 if (my_rank == 0){ 

    float*  f = new float[M*M];       // image
    float* R = new float [(M/p)*M];  // result
      cout<<"N= "<<N<<endl;
     for( i  = 0; i < M; i++){
        for( j = 1; j < M; j = j+2){
         *(f + i* M + j) = 0;
       }   
    }
   for( i  = 0; i < M; i++){
        for( j = 0; j < M; j = j + 2){
         *(f + i* M + j) = (W*W);
       }       
    }

   start = MPI_Wtime();   // starts monitoring execution time

 /* send partitions of the image for each spawned process */
   for(int k  = 1; k < p; k++){
    if(k == (p-1)){
    MPI_Send(f + ((M/p)*k-((W-1)/2)) * M, Nm * M, MPI_FLOAT, k, 0, MPI_COMM_WORLD);
    }
    else { MPI_Send(f + ((M/p)*k-((W-1)/2)) * M, N * M, MPI_FLOAT, k, 0, MPI_COMM_WORLD);}
     }
  
 /* rank 0 process convoluting */
  for( i = 0; i < (M/p); i++){
   for(  j = 0; j < M; j++){
     for (int k = 0; k < W; k++){
      for ( int l = 0; l < W; l++){
       int  row = i + (k - (W/2));
       int  col = j + (l - (W/2));
        if(row >= 0 && row < (M/p) && col >= 0 && col <  M){
            *(R + i* M + j) += *(f + row * M + col) * h[k][l];
            }
         }
      }
    }
  }

/* rank 0 process writes result into its partition of image matrix */
   for(int s = 0; s < (M/p); s++){
    for(int h = 0; h < M; h++){
     *(f + s * M + h) = *(R + s * M + h);
     }
   }

 delete [] R;

 /* rank 0 process receives and writes result into each partition */
  for (q = 1; q < p; q++){
  MPI_Recv(f + ((M/p)*q) * M,(M/p)*M, MPI_FLOAT, q, 0, MPI_COMM_WORLD, &status);
  }

 end = MPI_Wtime();
 cout <<"the parallel execution time for " << p << " process(es) is: "<< end - start <<"secs" <<endl;
 cout<<"  "<<endl;

     /* prints a section of the output image */
  cout<<"A section of the output image matrix: "<<endl;
  for(  i =1625; i < 1635; i++)
   {
    for(  j = 0; j < 9; j++)
    {

  cout<< *(f + i * 9 + j) <<" ";

    }
    cout<<endl;
   }
 cout<<" "<<endl;

delete [] f;
  
 } 
 
 else if(my_rank == (p-1)) {
   /* other  processes execute this block*/

         float* buff = new float [N*M];     // comm. buffer
         float* result = new float[(M/p)*M]; //result
  
      /* each receives its partition of the image from rank 0 process*/
    
     MPI_Recv(buff, Nm * M, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
         
        /* each process convolutes */
   for( i = 0; i < (M/p); i++){
   for(  j = 0; j < M; j++){
     for (int k = 0; k < W; k++){
      for ( int l = 0; l < W; l++){
       int  row = i + (k - (W/2));
       int  col = j + (l - (W/2));
        if(row >= ((W-3)/2) && row < (M/p) && col >= 0 && col <  M){
            *(result + i * M + j)  += *(buff + row * M + col) * h[k][l];
           }
        }
      }
    }
  }

 /* each processsends its result to rank 0 process */
 MPI_Send(result, (M/p)*M, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);

    /* each process deallocates memory */
delete [] buff;
delete [] result;
 }

else {
 float* buffn = new float [N*M];     // comm. buffer
 float* resultn = new float[(M/p)*M]; //result

 MPI_Recv(buffn, N * M, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);

     /* each process convolutes */
   for( i = 0; i < (M/p); i++){
   for(  j = 0; j < M; j++){
     for (int k = 0; k < W; k++){
      for ( int l = 0; l < W; l++){
       int  row = i + (k - (W/2));
       int  col = j + (l - (W/2));
        if(row >= ((W-3)/2) && row < (M/p) && col >= 0 && col <  M){
            *(resultn + i * M + j)  += *(buffn + row * M + col) * h[k][l];
           }
        }
      }
    }
  }

 /* each processsends its result to rank 0 process */
 MPI_Send(resultn, (M/p)*M, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);

    /* each process deallocates memory */
delete [] buffn;
delete [] resultn;

}

 MPI_Finalize();

return 0;
}
