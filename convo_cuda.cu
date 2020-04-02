#include <iostream>
#include <cstdlib>
#include <cuda.h>
#include <cmath>

#define M 2048
#define W 15
#define w 3
#define threshold 80

using namespace std;

 __global__  void smoothening_kernel(float* d_filter,float* d_raw_image,float* d_hx,float* d_hy,float* d_gx,float* d_gy,float* d_smooth_image,float* d_edged_image,int block_size){                  
               int Bx = blockIdx.x;
               int By = blockIdx.y;
               int Tx = threadIdx.x;
               int Ty = threadIdx.y;
     /* defining row and column index tp parse through filters and image*/
   int rowd = By* block_size + Ty;
   int columd = Bx* block_size + Tx;
     
     /*boundaries checking*/
   int rr = rowd - W/2;
   int cc = columd - W/2;
   float acc = 0.0;
    /*convolution for smmothening*/   
      for(int k = 0; k < W; k++ ){
       for(int l = 0; l < W; l++){
          if((rr + k) >= 0 && (rr + k) < M && (cc + l) >= 0 && (cc + l) < M){
           acc += d_raw_image[(rr + k) * M + (cc + l)] * d_filter[k * W + l];
        }  
      }
   d_smooth_image[rowd * M + columd] = acc;
    }

    /*convolution for edge detection */
     int mm = rowd - w/2;
     int nn = columd - w/2;
     float acc1 = 0.0;
     float acc2 = 0.0;
   for(int k = 0; k < w; k++ ){
    for(int l = 0; l < w; l++){
          if((mm + k) >= 0 && (mm + k) < M && (nn + l) >= 0 && (nn + l) < M){
           acc1 += d_smooth_image[(mm + k) * M + (nn + l)] * d_hx[k * w + l];
           acc2 += d_smooth_image[(mm + k) * M + (nn + l)] * d_hy[k * w + l];
         }
      } 
   d_gx[rowd * M + columd] = acc1;
   d_gy[rowd * M + columd] = acc2;
    }
 
  // gradient magnitude of spatial domains
    d_edged_image[rowd * M + columd] = sqrt(pow(d_gx[rowd * M + columd], 2) + pow(d_gy[rowd * M + columd], 2));
    if(d_edged_image[rowd * M + columd] > threshold){d_edged_image[rowd * M + columd]  = 255;}
      else{d_edged_image[rowd * M + columd] = 0;}
     }


 int main(int argc, char* argv[]){

int block_size = atoi(argv[1]);

float h_filter[W][W];  //Gaussian filter
float h_x[w][w] = {{-1.0,0.0,1.0},{-2.0,0.0,2.0},{-1.0,0.0,1.0}}; // Sobel operator
float h_y[w][w] = {{-1.0,-2.0,-1.0},{0.0,0.0,0.0},{1.0,2.0,1.0}}; //Sobel operator
 
 double sigma = 1.5;       
 float P = 1.0/(2* M_PI * sigma*sigma);
 float Q = 2.0* M_PI * sigma*sigma;
 float sum = 0.0;
 long image_size;
 size_t elements;
 int L = (W-1)/2;

  /*initializing gaussian filter*/
 for(int x = -W/2; x <= W/2; x++){
  for(int y = -W/2; y <= W/2; y++){
    int I = (x+ W/2) - L;
    int J = (y+ W/2) - L;
    h_filter[x + W/2][y + W/2] = P*(exp(-(I*I + J*J)/Q));
    sum += h_filter[x + W/2][y + W/2];
   }
  }

 for(int i = 0; i < W; i++){
  for(int j = 0; j < W; j++){
  h_filter[i][j]/= sum;
  }
 }
  // verify gaussian filter
  cout<<"guassian filter" <<endl;
  for(int q = 0; q < 15; q++){
   for(int z = 0; z <15; z++){
     cout<<h_filter[q][z]<<" ";
    }
 cout<<endl;
  }
 cout<<" "<<endl;

FILE* fp_in, *fp_out1, *fp_out2;
   fp_in = fopen ("Rainier2048_noise.bin","rb");
      if(fp_in == NULL){cout<<"FILE ERROR!"<<endl;
          exit(1); }

   //obtain file size
  fseek(fp_in, 0, SEEK_END);
  image_size = ftell(fp_in);
  rewind(fp_in);

  // allocate buffer of image size
  unsigned char* buffer = (unsigned char*)malloc(sizeof(unsigned char) * image_size);
  unsigned char* buffer1 = (unsigned char*)malloc(sizeof(unsigned char) * image_size);
 //copy file into buffer
 elements = fread(buffer, sizeof(unsigned char), image_size, fp_in);
 if(elements != image_size){cout<<"READ ERROR! "<<endl;
       exit(2);}

  fclose(fp_in);

  float* fptr = (float*)malloc(sizeof(float)* M * M);
  
  //typecast from char to float
   for(int row = 0; row < M; row++){
   for(int col = 0; col < M; col++){
    fptr[row * M + col] = (float) buffer[row * M + col];
  }
 }

   cout<<"raw image" <<endl;
  for(int q = 1024; q < 1034; q++){
   for(int z = 1525; z <1535; z++){
     cout<<buffer[q * M + z]<<" ";
    }
 cout<<endl;
  }

  cout<<"raw image of float type" <<endl;
  for(int q = 1024; q < 1034; q++){
   for(int z = 1525; z <1535; z++){
     cout<<fptr[q * M + z]<<" ";
    }
 cout<<endl;
  }
cout<<" "<<endl;
 
 float* smooth_image = (float*)malloc(sizeof(float)* M * M);
 float* edged_image = (float*)malloc(sizeof(float)* M * M); 

   
  float* d_gx;
  float* d_gy;
  float* d_hx;
  float* d_hy;
  float* d_raw_image;
  float* d_filter;
  float* d_smooth_image;
  float* d_edged_image;

  cudaMalloc((void**)&d_hx,sizeof(float)* w * w);
  cudaMalloc((void**)&d_hy,sizeof(float)* w * w);
  cudaMalloc((void**)&d_filter,sizeof(float)* W * W);
  cudaMalloc((void**)&d_raw_image,sizeof(float)* M * M);
  cudaMalloc((void**)&d_smooth_image,sizeof(float)* M * M);
  cudaMalloc((void**)&d_edged_image,sizeof(float)* M * M);
  cudaMalloc((void**)&d_gx,sizeof(float)* M * M);
  cudaMalloc((void**)&d_gy,sizeof(float)* M * M);
 
 /* measuring execution time */
 cudaEvent_t start, stop;
 cudaEventCreate(&start);
 cudaEventCreate(&stop);

  /*copy image and filters from host to device */
   cudaMemcpy(d_raw_image, fptr, sizeof(float) * M * M, cudaMemcpyHostToDevice);
   cudaMemcpy(d_filter,h_filter , sizeof(float) * W * W, cudaMemcpyHostToDevice);
   cudaMemcpy(d_hx, h_x , sizeof(float) * w * w, cudaMemcpyHostToDevice);
   cudaMemcpy(d_hy, h_y , sizeof(float) * w * w, cudaMemcpyHostToDevice);

  /*define block size and grid size and invoke kernel*/
    dim3 threadsPerBlock(block_size, block_size);
    int numblocks = M / block_size;
    dim3 blocksPerGrid(numblocks, numblocks);

   cudaEventRecord(start);
 smoothening_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_filter,d_raw_image,d_hx,d_hy,d_gx,d_gy,d_smooth_image,d_edged_image,block_size);

     /* copy results from device to host */
  cudaMemcpy(smooth_image, d_smooth_image, sizeof(float) * M * M, cudaMemcpyDeviceToHost);
  cudaMemcpy(edged_image, d_edged_image, sizeof(float) * M * M, cudaMemcpyDeviceToHost);

 cudaEventRecord(stop);
 cudaEventSynchronize(stop);

float milliseconds = 0.0;

 cudaEventElapsedTime(&milliseconds, start, stop);

 cout<<"he parallel execution time for block size "<< block_size << " is "<< milliseconds <<" secs" << endl;

  /* write edge detected image to file*/
   for(int row = 0; row < M; row++){
   for(int col = 0; col < M; col++){
    buffer[row * M + col] = (unsigned char) smooth_image[row * M + col];
    buffer1[row * M + col] = (unsigned char) edged_image[row * M + col];
   }
  }
     cout<<"smoothened_image buffered"<<endl;
    for(int ir = 1024; ir < 1034; ir++){
      for(int ic = 1525; ic < 1535; ic++){
    cout<< *(buffer + ir * M + ic) <<" ";
   }
 cout<<endl;
 }
 cout<<" "<<endl;

 
  
   fp_out1 = fopen("smoothened_image_cuda.bin", "wb");
   fwrite(buffer, sizeof(unsigned char), image_size, fp_out1);
    fclose(fp_out1);

   fp_out2 = fopen("Edge_detected_image_cuda.bin", "wb");
   fwrite(buffer1,sizeof(unsigned char), image_size, fp_out2);
    fclose(fp_out2);
 
   cout<<"smoothened image" <<endl;
  for(int q = 1024; q < 1034; q++){
   for(int z = 1525; z <1535; z++){
     cout<<smooth_image[q * M + z]<<" ";
    }
 cout<<endl;
  }
cout<<" "<<endl;

   cout<<"edged_image buffered"<<endl;
  for(int ir = 1024; ir < 1034; ir++){
    for(int ic = 1525; ic < 1535; ic++){
   cout<< *(buffer1 + ir * M + ic) <<" ";
   }
 cout<<endl;
 }
 cout<<" "<<endl;

   cout<<"edged_image" <<endl;
  for(int q = 1024; q < 1034; q++){
   for(int z = 1525; z <1535; z++){
     cout<<edged_image[q * M + z]<<" ";
    }
 cout<<endl;
  }
cout<<" "<<endl;

/* free device memory*/
cudaFree(d_raw_image);
cudaFree(d_hx);
cudaFree(d_hy);
cudaFree(d_smooth_image);
cudaFree(d_edged_image);
cudaFree(d_gx);
cudaFree(d_gy);
cudaFree(d_filter);

/*free host memory*/
delete[] fptr;
delete[] smooth_image;
delete[] buffer;
delete[] buffer1;
delete[] edged_image;


 return 0;
 }
