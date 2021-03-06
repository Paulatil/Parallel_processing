#include <iostream>
#include <cstdlib>
#include <cuda.h>
#include <ctime>
#include <sys/time.h>

using namespace std;

     /*the kernel code to run on the GPU device */
    __global__ 
          void matrix_mult_kernel(float* A, float* B, float* C, int M, int block_size){
           /* using specified conventions*/
               int Bx = blockIdx.x;
               int By = blockIdx.y;
               int Tx = threadIdx.x;
               int Ty = threadIdx.y;
               int Grid_size = M / block_size;
             
             
      /* defining row and column index tp parse through matrix A & B */
   int rowd = (By * block_size) + Ty;
   int columd = (Bx * block_size) + Tx;
   int rowds;
   int columds;
   extern  __shared__ float smem[];
                       float* ads = (float*)(smem);
                       float*  bds =(float*)(smem + block_size * block_size) ;
                       
  float tempsum = 0;

        ads[Ty * block_size + Tx] = 0;
        bds[Ty * block_size + Tx] = 0;
        

  

        /* copy data into the allocated shared memory */
        for(int k = 0; k < Grid_size; k++){
         rowds = rowd * M + (k * block_size + Tx);
         columds = (k * block_size + Ty) * M + columd;
        ads[Ty * block_size + Tx] =  A[rowds];
        bds[Ty * block_size + Tx] = B[columds];   
    
     __syncthreads();   // synchronize 

     for(int i = 0; i < block_size; i++)
     tempsum += ads[Ty * block_size + i] * bds[i * block_size + Tx];    
     
   __syncthreads();         
     }
   C[rowd * M + columd] = tempsum;
   
  } 
  


 int main(int argc, char* argv[]){ 
  int M = 4096;
  int B = atoi(argv[1]);   //block size
      
 
 /*allocate matrixes A, B, C in host memory*/ 
  float* ahptr = (float*)malloc(sizeof(float)* M * M);
  float* bhptr = (float*)malloc(sizeof(float)* M * M);
  float* chptr = (float*)malloc(sizeof(float)* M * M);
  float* dhptr = (float*)malloc(sizeof(float)* M * M);
 /* initialize matrices a, b in host memory*/
  
  for(int i = 0; i < M; i++){
   for(int j = 0; j < M; j++){
     *(ahptr + i * M + j) = ((i+1)*(j+1))/(float)M;
     *(bhptr + i * M + j) = (float)(j+1)/(i+1); 
     *(chptr +i * M + j) = 0;
     *(dhptr +i * M + j) = (i+1)*(j+1);
   }
 }
   
     //verify result   
   cout<<"result verifier"<<endl;
    for(int w = 2044; w < 2052; w++){
     for(int s = 0; s < 8; s++){
     cout<< *(dhptr + w * M + s)<<" ";
     }
   cout<<endl;
   }
  cout<<" "<<endl;

    /*allocate memoryon the device*/
 
   float* ad;
   float* bd;
   float* cd;

  cudaMalloc((void**)&ad,sizeof(float)* M * M);
  cudaMalloc((void**)&bd,sizeof(float)* M * M); 
  cudaMalloc((void**)&cd,sizeof(float)* M* M);
   /*set shared memory to 48kB and Li cache to 16kB*/
  cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);  
   /* measuring the execution time*/
 cudaEvent_t start, stop;
 cudaEventCreate(&start);
 cudaEventCreate(&stop);
    /*copy matrices from host to device */
   cudaMemcpy(ad, ahptr, sizeof(float) * M * M, cudaMemcpyHostToDevice);
   cudaMemcpy(bd, bhptr, sizeof(float) * M * M, cudaMemcpyHostToDevice);

     /*invoking the kernel */
     int block_size = B;
    dim3 threadsPerBlock(block_size, block_size);
    int numblocks = M / block_size;
    dim3 blocksPerGrid(numblocks, numblocks);

 cudaEventRecord(start);
 matrix_mult_kernel<<< blocksPerGrid, threadsPerBlock, sizeof(float) * block_size * block_size * 2 >>>(ad, bd, cd, M, B);

   /* copy result from device to host */
 cudaMemcpy(chptr, cd, sizeof(float) * M * M, cudaMemcpyDeviceToHost);

 cudaEventRecord(stop);
 cudaEventSynchronize(stop);

float milliseconds = 0.0;
 
 cudaEventElapsedTime(&milliseconds, start, stop);

 cout<<"the parallel execution time for block size "<< B << " is "<< milliseconds <<" msecs" <<endl; 

    /*print a section of the result to verify result*/ 
   cout<<" a section of the GPU result;"<<endl;
    for(int h = 2044; h < 2052; h++){
   for(int t = 0; t < 8; t++){
   cout<< *(chptr + h * M + t) <<" ";
  }
  cout<<endl;
 }
 cout<<" "<<endl;

 /* free device memory*/
cudaFree(ad);
cudaFree(bd);
cudaFree(cd);

 /*free host memory*/
delete [] ahptr;
delete [] bhptr;
delete [] chptr;

 return 0;
 
 };
