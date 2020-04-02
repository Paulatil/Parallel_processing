#include <iostream>
#include <cstdio>
#include <omp.h>
#include <cmath>
#include <cstdlib>
#include <ctime>
#define M 3240

using namespace std;

/**************** GLOBAL MATRICES A< B< C ****************************/
 float A[M][M];
 float B[M][M];
 float C[M][M];
 float K, error;

/////////////////////////////////////////////////////
 double omp_get_wtime();
//-----------------------------------------------------------//
void omp_set_num_threads(int nthreads);
//////////// matrix multip;ication function /////////////////////////////
 void matrix_multiply(float A[][M], float B[][M], float C[][M]){
  float result;
 int i, j, k; 
 #pragma omp parallel
    { 
    #pragma omp for  schedule(static)  private(i,j,k)  reduction(+: result) collapse(2)
     for( i = 0; i < M; i++){
      for( j = 0; j < M; j++){
       for( k = 0; k < M; k++ ){
         result+= A[i][k] * B[k][j];
         }
     C[i][j] = result;
     }
   }
  cout<<"the number of threads in matrix_multiplication: "<<omp_get_num_threads()<<endl;
 }
}
/////////////// AVERAGE PERCENT ERROR /////////////////////////////////
 void Average_Percent_Error(float C[][M]){
  int i, j;
  #pragma omp parallel  
   {
   #pragma omp for  schedule(static) private(i, j) reduction(+: K) collapse(2)
  for ( i = 0; i < M; i++){
    for( j = 0; j < M; j++){
         K += (C[i][j] - (i+1)*(j+1))/((i+1)*(j+1));
      }
    }
    cout<<"the number of threads in error compuation: "<<omp_get_num_threads()<< endl; 
  
        
  #pragma omp master 
       {
     error = K *(100/(M*M));
      }
    }
   cout <<"the K value is: "<< K << endl;
   cout<<"the AVE is: "<< error << endl;  
}

/////////////// MAIN FUNCTION ////////////////////////////////////
 int main(){
 double start_time;
double end_time;
int i, j;
int nthreads;
   
   start_time = omp_get_wtime();
    
    #pragma omp parallel 
    {
    #pragma omp for schedule(static) private(i,j) collapse(2)
   for (i = 0; i < M; i++){
    for(j = 0; j < M; j++){
     A[i][j] = ((i+1)*(j+1)/M);
     B[i][j] = (j+1)/(i+1);
     C[i][j] = 0;
      }
    }
  cout<<"the number of threads for matrix A, B, C initializations: "<< omp_get_num_threads() <<endl;  
  nthreads = omp_get_num_threads(); 
 } 

 matrix_multiply(A, B, C);   
 Average_Percent_Error(C); 
 
 end_time = omp_get_wtime(); 
  cout<<"the execution time for "<<nthreads<<" is:"<< end_time - start_time << endl;

  cout<<" A 6*6 centre block of matrix C is:"<<endl;

 for(i = 1616; i < 1623; i++){
   for( j = 1616; j < 1623; j++){
    cout<< C[i][j] << " ";
  }
  cout<< endl;
 }


return 0;
}  
