#include <iostream>
#include <cmath>
#define M  4
using namespace std;

/*************** declaring matrices A, B, C as global 2D arrays********/
float A[M][M];
float B[M][M];
float C[M][M];

///////////////////////////////////////////////////////////////////////
/* matrix-multiplication****************************************/
void  matrix_mult(float matrixA[][M], float matrixB[][M],float matrixC[][M])
{
  for(int i = 0; i < M; i++)
    {
    for(int j = 0; j < M; j++)
     {
       for(int k  = 0; k < M; k++)
         {
        C[i][j] += matrixA[i][k] * matrixB[k][j];

          }
       }
     }
  cout<<endl;
  }
/////////////////////////////////////////////////////////////////////
/**********printing each row of the matrices A,B,C ****************/
void print_row_by_row(float matrixA[][M], float matrixB[][M], float matrixC[][M])
  {
   for(int i = 0; i < M; ++i)
    { cout<<"elements of row "<< i + 1 <<" of matrix C are: ";
     for(int j = 0; j < M; ++j)
      {
        cout<< C[i][j]<<" ";
          if(j == M-1){
             cout<<""<<endl;
               }
       }
 cout<<"--------------------------------------------"<<endl;

   for(int i = 0; i < M; ++i)
    { cout<<"elements of row "<< i + 1 <<" of matrix A are: ";
     for(int j = 0; j < M; ++j)
      {
        cout<< A[i][j]<<" ";
          if(j == M-1){
             cout<<""<<endl;
               }
       }
     }
    cout<<"-----------------------------------------------"<<endl;

 for(int i = 0; i < M; ++i)
    { cout<<"elements of row "<< i + 1 <<" of matrix B are: ";
     for(int j = 0; j < M; ++j)
      {
        cout<< B[i][j]<<" ";
          if(j == M-1){
             cout<<""<<endl;
               }
       }
     }
    cout<<"---------------------------------------------"<<endl;

 }
//////////////////////////////////////////////////////////////////////////
/***** teh main function ***********************/
int main()
{

/*****Initializing matrices A,B,C***************/
   for(int i = 0; i < M; i++){
     for(int j = 0; j < M; j++){
       C[i][j] = 0;
       }
     }
////////////////////////////////////////////////////////
/*getting and storing elements into matrices A & B */
///////////////////////////////////////////////////////
   cout<<"Enter the elements of matrix A:"<< endl;
    for(int row = 0; row < M; row++){
     for(int column = 0; column < M; column++){
        cout<<"enter element A"<< row + 1  << column + 1  << ":";
        cin >> A[row][column];
       }
     }
    cout<<endl;

   cout<<"Enter the elements of matrix B:"<<endl;
    for(int row = 0; row < M; row++ ){
      for(int column = 0; column < M; column++){
       cout<<"enter element B"<< row + 1 << column + 1 <<":";
       cin >> B[row][column];
       }
     }
    cout<<endl;


     matrix_mult(A,B,C);
  print_row_by_row(A,B,C);

return 0;

 }

