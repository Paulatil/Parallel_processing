#include <iostream>
#include <cmath>
#include <omp.h>
#define M 8
#define W 3

using namespace std;

/******global Image, filter and filtered Image matrices*************/
float f[M][M]; //Image
float g[M][M]; //filter
float h[W][W]; //filtered image

/********2D convolution function********************/
void convolute(float image[][M],float filter[][W],float Fimage[][M])
{
 #pragma omp parallel for 
   
 for(int i = 0; i < M; i++)
 {
   for(int j = 0; j < M; j++) 
   {
     for(int k = 0; k < W; k++ )
     {
       for(int l = 0; l < W; l++)
        {
          int row = i + (k - (W/2));  // new index for boundary checking of image matrix
          int col = j + (l - (W/2));  
        if(row >= 0 && row < M && col >= 0 && col < M){
            g[i][j] += f[row][col] * h[k][l];} 
         }                                      
      }
    }   
  }

 }



/*********print function for outputting elements in each matrix********************/
void print_matrix(float image[][M],float filter[][W],float Fimage[][M])
{ 
 #pragma omp master
 {
   cout<<"the filtered image matrix:"<<endl;
  for(int i = 0; i < M; i++)
   {
    for(int j = 0; j < M; j++)
    {
     cout<< Fimage[i][j] <<" ";
     
    }
    cout<<endl;
   }
   cout<<"-----------------------------------"<<endl;
///////////////////////////////////////////////
 cout<<"the input image matrix:"<<endl;
  for(int i = 0; i < M; i++)
   {
    for(int j = 0; j < M; j++)
    { 
     cout<< image[i][j] <<" ";
    }
    cout<<endl;
   }
   cout<<"-----------------------------------"<<endl;
////////////////////////////////////////////////
 cout<<"the filter matrix:"<<endl;
  for(int i = 0; i < W; i++)
   {
    for(int j = 0; j < W; j++)
    {
     cout<< filter[i][j] <<" ";
    }
   cout<<endl; 
   }
   cout<<endl;
 
 }
}
/****************** amin function *************************/
int main()
{
 #pragma omp parallel for
 
  for(int i = 0; i < M; i++)
  {
   for(int j = 0; j < M; j++)
    {
       g[i][j] = 0;
       f[i][j] = 0;        
    } 
  }


 #pragma omp parallel for
 
    for(int row = 0; row < M; row++)
    {
      for(int column = 0; column < M; column+=2)
    {
       f[row][column] = 841;  
      }
     }
   

 #pragma omp parallel for 
 
   
    for(int frow = 0; frow < W; frow++)
   {
    for(int fcol = 0; fcol < W; fcol++)
    {
     h[frow][fcol] = 1.0/9;

     }
   }
 
   convolute(f, h, g);
   print_matrix(f, h, g);


return 0;
}
