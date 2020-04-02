#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#define W 29
#define M 336

using namespace std;
//////////// global vabriables ////////////////
float f[M][M]; //Image and filtered image
float h[W][W]; //filter

/////////////// main function ///////
int main(){
int G = 336;
double start_time, end_time;
int nthreads;


 start_time = omp_get_wtime();

#pragma omp parallel for
 for ( int i = 0; i < M; i++){
   for ( int j = 0; j < M; j++){
    f[i][j] = 0;
   }
 }

#pragma omp parallel for
 for ( int i = 0; i < M; i++){
   for ( int j = 0; j < M; j+=2){
    f[i][j] = 841;
  }
 }


#pragma omp parallel for
 for ( int i = 0; i < W; i++){
   for (int j = 0; j < W; j++){
     h[i][j] = 1.0/(W*W);
  }
 }

////////// convolution done in this block of the code  ////////////////////
#pragma omp parallel
 {
 nthreads = omp_get_num_threads();  
   int tid = omp_get_thread_num();
   #pragma omp for
    for(int start = 0; start < nthreads; start++){
   float** temp = new float* [(G/nthreads)];
   for(int s = 0; s < (G/nthreads); s++){
    temp[s] = new float [G];}
 
  for(int i = (tid*(G/nthreads)); i < ((tid*(G/nthreads))+(G/nthreads)); i++){
   for( int j = 0; j < G; j++){
     for (int k = 0; k < W; k++){
      for ( int l = 0; l < W; l++){
      int row = i + (k - (W/2));
       int col = j + (l - (W/2));
        if(row >= (tid*(G/nthreads)) && row < ((tid*(G/nthreads))+(G/nthreads)) && col >= 0 && col <  G){ 
            temp[i][j] += f[row][col] * h[k][l];
               }
              }
             }
           }
         }
       for(int i = (tid*(G/nthreads)); i < ((tid*(G/nthreads))+(G/nthreads)); i++){
   for( int j = 0; j < G; j++){
       f[i][j] = temp[i][j];
         }
       }  
    for(int t = 0; t < (G/nthreads); t++){
    delete[] temp[t];
    delete[] temp;
    }
   }
  }
 cout<<"Number of threads N = "<< nthreads <<endl;

 end_time = omp_get_wtime();
cout<<"the parallel execution time is "<< end_time - start_time << endl;
//////////////print out the output image ///////////////
 #pragma omp master 
{
  cout<<"the filtered image matrix:"<<endl;
  for( int i = 163; i < 172; i++)
   {
    for(int j = 0; j < 17; j++)
    {

  cout<< f[i][j] <<" ";

    }
    cout<<endl;
   }
 } 
 

return 0;

}


