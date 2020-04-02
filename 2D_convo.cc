#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#define M 3360
#define W 29

using namespace std;
//////////// global vabriables ////////////////

float f[M][M]; //Image
float h[W][W]; //filter
float g[M][M]; //filtered image

/////////////// main function ///////
int main(){

double start_time, end_time;
int nthreads;


 start_time = omp_get_wtime();

#pragma omp parallel for
 for ( int i = 0; i < M; i++){
   for ( int j = 0; j < M; j++){
    f[i][j] = 0;
    g[i][j] = 0;
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
#pragma omp parallel for
 for(int i = 0; i < M; i++){
   for( int j = 0; j < M; j++){
     for (int k = 0; k < W; k++){
      for ( int l = 0; l < W; l++){
      int row = i + (k - (W/2));
       int col = j + (l - (W/2));
        if(row >= 0 && row < M && col >= 0 && col <  M){ 
          g[i][j] += f[row][col] * h[k][l];
          }
        }
      }
    }
nthreads = omp_get_num_threads();
}

 cout<<"Number of threads N = "<< nthreads <<endl;

 end_time = omp_get_wtime();
cout<<"the parallel execution time is "<< end_time - start_time << endl;
//////////////print out the output image ///////////////
 #pragma omp master 
{
  cout<<"the filtered image matrix:"<<endl;
  for( int i = 1675; i < 1685; i++)
   {
    for(int j = 0; j < 17; j++)
    {

  cout<< g[i][j] <<" ";

    }
    cout<<endl;
   }
 }
 
 #pragma omp single
  {
     cout<<"the image matrix:"<<endl;
  for( int i = 1675; i < 1685; i++)
   {
    for( int j = 0; j < 11; j++)
    {

  cout<< f[i][j] <<" ";

    }
    cout<<endl;
   }
 }

#pragma omp single 
{
cout<<"the filter matrix:"<<endl;
  for( int i = 0; i < 11; i++)
   {
    for( int j = 0; j < 11; j++)
    {

  cout<< h[i][j] <<" ";

    }
    cout<<endl;
   }
}

return 0;

}


