#include <iostream>
#include <cmath>
#include <cstdlib>
#include <mpi.h>

#define W 15
#define w 3
#define M 2048
#define threshold 80
using namespace std;

/* GLOBAL VARIABLES */
float h_filter[W][W];  //Gaussian filter
float h_x[w][w] = {{-1.0,0.0,1.0},{-2.0,0.0,2.0},{-1.0,0.0,1.0}}; // Sobel operator
float h_y[w][w] = {{-1.0,-2.0,-1.0},{0.0,0.0,0.0},{1.0,2.0,1.0}}; //Sobel operator

/* MAIN FUNCTION*/
int main(int argc, char* argv[]){


 double start, end, sigma = 1.5;             /* timestamps */
 int i, j, q, I, J;
 int L = (W-1)/2;
 int ls = (w-1)/2;
 float P = 1.0/(2* M_PI * sigma*sigma);
 float Q = 2.0* M_PI * sigma*sigma;
 float sum = 0.0;
 MPI_Status status;

//generate  gaussian filter
 for(int x = -W/2; x <= W/2; x++){
  for(int y = -W/2; y <= W/2; y++){
     I = (x+ W/2) - L;
     J = (y+ W/2) - L;
    h_filter[x + W/2][y + W/2] = P*(exp(-(I*I + J*J)/Q));
    sum += h_filter[x + W/2][y + W/2];
   }
}
 
 for(i = 0; i < W; i++){
  for(j = 0; j < W; j++){
  h_filter[i][j]/= sum;  
  }
}

 MPI_Init(&argc, &argv);

 int my_rank, p;

 MPI_Comm_size(MPI_COMM_WORLD, &p);

 MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

 int N = (M/p) + (2*L);       /* row size of comm buffer */
 int Nm = (M/p) + L;  /* size of partition for last process */
 int Ns = (M/p) + (2*ls); // partition for edge detection
 int Nms = (M/p) + ls;    // parition for edge detection
/* rank 0 process generates image matrix f */
 if (my_rank == 0){
    float* f = (float*)malloc(sizeof(float)* M * M);//allocating memory fo image
    float* fini = (float*)malloc(sizeof(float)* M * M);
    float* R = new float[(M/p) * M];  // result of filtering
    float* g_x = new float[(M/p) * M];
    float* g_y = new float[(M/p) * M];
    unsigned char* buffer;
   unsigned  char* buffer1;
   long image_size;
   size_t elements; 
    /*read data from binary file*/
     FILE* fp_in, *fp_out1, *fp_out2;
      fp_in = fopen ("Rainier2048_noise.bin","rb");
      if(fp_in == NULL){cout<<"FILE ERROR!"<<endl;
          exit(1); }
   
   //obtain file size
  fseek(fp_in, 0, SEEK_END);
  image_size = ftell(fp_in);
  rewind(fp_in);
 
  // allocate buffer of image size
  buffer = (unsigned char*)malloc(sizeof(unsigned char) * image_size);
  buffer1 = (unsigned char*)malloc(sizeof(unsigned char) * image_size);
 //copy file into buffer
 elements = fread(buffer, sizeof(unsigned char), image_size, fp_in);
 if(elements != image_size){cout<<"READ ERROR! "<<endl;
       exit(2);}
 
  fclose(fp_in);

// typecast from char to float
 for(int row = 0; row < M; row++){
   for(int col = 0; col < M; col++){
    f[row * M + col] = (float) buffer[row * M + col];
  }
 }
  //print a section of buffer
  for(int ir = 2030; ir < M; ir++){
  for(int ic = 2005; ic < 2023; ic++){
   cout<< *(buffer + ir * M + ic) <<" "; 
   }
 cout<<endl;
 }
 cout<<" "<<endl;

  start = MPI_Wtime();   // starts monitoring execution time

  /* send partitions of raw image for each spawned process */
   for(int k  = 1; k < p; k++){
    if(k == (p-1)){
    MPI_Send(f + ((M/p) * k - L) * M, Nm * M, MPI_FLOAT, k, 0, MPI_COMM_WORLD);
    }
    else { MPI_Send(f + ((M/p) * k - L) * M, N * M, MPI_FLOAT, k, 0, MPI_COMM_WORLD);}
     }

  /* rank 0 process convoluting */
  for( i = 0; i < (M/p); i++){
   for(  j = 0; j < M; j++){
     for (int k = 0; k < W; k++){
      for ( int l = 0; l < W; l++){
       int  row = i + (k - (W/2));
       int  col = j + (l - (W/2));
        if(row >= 0 && row < (M/p) && col >= 0 && col <  M){
            *(R + i* M + j) += *(f + row * M + col) * h_filter[k][l];}
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

/* rank 0 process collates results for smoothening operation */
 for (q = 1; q < p; q++){
  MPI_Recv(f + ((M/p)*q) * M,(M/p)*M, MPI_FLOAT, q, 0, MPI_COMM_WORLD, &status);
  }
  
 /* send partitions of smoothened image for edge detection operation */
   for(int k  = 1; k < p; k++){
    if(k == (p-1)){
    MPI_Send(f + ((M/p) * k - ls) * M, Nms * M, MPI_FLOAT, k, 0, MPI_COMM_WORLD);
    }
    else { MPI_Send(f + ((M/p) * k - ls) * M, Ns * M, MPI_FLOAT, k, 0, MPI_COMM_WORLD);}
     }

 /*  edge detection process */
  for( i = 0; i < (M/p); i++){
   for(  j = 0; j < M; j++){
     for (int k = 0; k < w; k++){
      for ( int l = 0; l < w; l++){
       int  row = i + (k - (w/2));
       int  col = j + (l - (w/2));
        if(row >= 0 && row < (M/p) && col >= 0 && col <  M){
            *(g_x + i* M + j) += *(f + row * M + col) * h_x[k][l];
            *(g_y + i* M + j) += *(f + row * M + col) * h_y[k][l];
               }
             }
           }
         }
       }
   //gradient magnitude of spatial domains 
  for( i = 0; i < (M/p); i++){
   for(  j = 0; j < M; j++){
    *(R + i* M + j) = sqrt(pow(*(g_x + i* M + j), 2) + pow(*(g_y + i* M + j), 2)); 
    }
  }

 // thresholding gradient magnitude to preserve image
 for(int row = 0; row < (M/p); row++){
   for( int col = 0; col < M; col++){
     if(*(R + row * M + col) > threshold){*(R + row * M + col) = 255;}
      else {*(R + row * M + col) = 0;}
   }
  }
 /* rank 0 process writes result into its partition of image matrix */
   for(int s = 0; s < (M/p); s++){
    for(int h = 0; h < M; h++){
     *(fini + s * M + h) = *(R + s * M + h);
     }
   }
 /* rank 0 process receives final results from edge detection operation */
 for (q = 1; q < p; q++){
  MPI_Recv(fini + ((M/p)*q) * M,(M/p)*M, MPI_FLOAT, q, 0, MPI_COMM_WORLD, &status);
  }


 end = MPI_Wtime();
 cout <<"the parallel execution time for " << p << " process(es) is: "<< end - start <<" secs" <<endl;
 cout<<"  "<<endl;
  
  cout<<"A section of  smoothened_image: "<<endl;
  for(  i =1024; i < 1034; i++)
   {
    for(  j = 1525; j < 1535; j++)
    {

  cout<< *(f + i * M + j) <<" ";
    }
    cout<<endl;
   }
 cout<<" "<<endl;

/* write smoothened image to file*/
   for(int row = 0; row < M; row++){
   for(int col = 0; col < M; col++){
    buffer[row * M + col] = (unsigned char) f[row * M + col];
  }
 }

   fp_out1 = fopen("smoothened_image2048.bin", "wb");
   fwrite(buffer, sizeof(unsigned char), image_size, fp_out1);
    fclose(fp_out1);

    for(int ir = 1024; ir < 1034; ir++){
     for(int ic = 1525; ic < 1535; ic++){
   cout<< *(buffer + ir * M + ic) <<" ";
   }
 cout<<endl;
 }
 cout<<" "<<endl;

 
 cout<<"A section of edge_detected_image: "<<endl;
  for(  i =1024; i < 1034; i++)
   {
    for(  j = 1525; j < 1535; j++)
    {

  cout<< *(fini + i * M + j) <<" ";
    }
    cout<<endl;
   }
 cout<<" "<<endl;

  /* write edge detected image to file*/
   for(int row = 0; row < M; row++){
   for(int col = 0; col < M; col++){
    buffer1[row * M + col] = (unsigned char) fini[row * M + col];
  }
 }

   fp_out2 = fopen("Edge_detected_image2048.bin", "wb");
   fwrite(buffer1,sizeof(unsigned char), image_size, fp_out2);
    fclose(fp_out2);

   for(int ir = 1024; ir < 1034; ir++){
    for(int ic = 1525; ic < 1535; ic++){
   cout<< *(buffer1 + ir * M + ic) <<" ";
   }
 cout<<endl;
 }
 cout<<" "<<endl;

delete[] f;
delete[] fini;
delete[] R;
delete[] buffer;
delete[] buffer1;
delete[] g_x;
delete[] g_y;
 }

 else if(my_rank == (p-1)) {
   /* other  processes execute this block*/

         float* buff = new float[Nm*M];     // comm. buffer
         float* result = new float[(M/p)*M]; //result
         float* g_x = new float[(M/p) * M];
         float* g_y = new float[(M/p) * M];
       
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
            *(result + i * M + j)  += *(buff + row * M + col) * h_filter[k][l];
           }
        }
      }
    }
  }

 /* each processsends its result of smoothening to rank 0 process */
 MPI_Send(result, (M/p)*M, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);

   float* buff1 = new float[Nms*M]; //allocate memory

  /* each receives its partition of smoothened image for edge detection operation */
     MPI_Recv(buff1, Nms * M, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);


 /*  edge detection process */
  for( i = 0; i < (M/p); i++){
   for(  j = 0; j < M; j++){
     for (int k = 0; k < w; k++){
      for ( int l = 0; l < w; l++){
       int  row = i + (k - (w/2));
       int  col = j + (l - (w/2));
        if(row >= 0 && row < (M/p) && col >= 0 && col <  M){
            *(g_x + i* M + j) += *(buff1 + row * M + col) * h_x[k][l];
            *(g_y + i* M + j) += *(buff1 + row * M + col) * h_y[k][l];
               }
             }
           }
         }
       }
  
  //gradient magnitude of spatial domains 
    for( i = 0; i < (M/p); i++){
      for(  j = 0; j < M; j++){
        *(result + i* M + j) = sqrt(pow(*(g_x + i* M + j), 2) + pow(*(g_y + i* M + j), 2));
         }
      } 
    // thresholding gradient magnitude to preserve image
      for(int row = 0; row < (M/p); row++){
       for( int col = 0; col < M; col++){
        if(*(result + row * M + col) > threshold){*(result + row * M + col) = 255;}
          else {*(result + row * M + col) = 0;}
          }
        }
                                       
 /* each processsends its result to rank 0 process */
 MPI_Send(result, (M/p)*M, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);

    /* each process deallocates memory */
delete [] buff;
delete[]buff1;
delete [] result;
delete[] g_x;
delete[] g_y;

 }

else {
 float* buffn = new float[N*M];     // comm. buffer
 float* resultn = new float[(M/p)*M]; //result
 float* g_x = new float[(M/p) * M];
 float* g_y = new float[(M/p) * M];

 MPI_Recv(buffn, N * M, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);

     /* each process convolutes */
   for( i = 0; i < (M/p); i++){
   for(  j = 0; j < M; j++){
     for (int k = 0; k < W; k++){
      for ( int l = 0; l < W; l++){
       int  row = i + (k - (W/2));
       int  col = j + (l - (W/2));
        if(row >= ((W-3)/2) && row < (M/p) && col >= 0 && col <  M){
            *(resultn + i * M + j)  += *(buffn + row * M + col) * h_filter[k][l];
           }
        }
      }
    }
  }

 /* each processsends its result of image smmothening to rank 0 process */
 MPI_Send(resultn, (M/p)*M, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);

 float* buffn1 = new float[Ns*M];  //allocate memory

 /*each process recieves a paritition of smoothened image for edge detection operation */
  MPI_Recv(buffn1, Ns * M, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
 /*  edge detection process */
  for( i = 0; i < (M/p); i++){
   for(  j = 0; j < M; j++){
     for (int k = 0; k < w; k++){
      for ( int l = 0; l < w; l++){
       int  row = i + (k - (w/2));
       int  col = j + (l - (w/2));
        if(row >= 0 && row < (M/p) && col >= 0 && col <  M){
            *(g_x + i* M + j) += *(buffn1 + row * M + col) * h_x[k][l];
            *(g_y + i* M + j) += *(buffn1 + row * M + col) * h_y[k][l];
               }
             }
           }
         }
       }

 //gradient magnitude of spatial domains 
  for( i = 0; i < (M/p); i++){
   for(  j = 0; j < M; j++){
   *(resultn + i* M + j) = sqrt(pow(*(g_x + i* M + j), 2) + pow(*(g_y + i* M + j), 2));
      }
    }
 
   // thresholding gradient magnitude to preserve image
   for(int row = 0; row < (M/p); row++){
    for( int col = 0; col < M; col++){
     if(*(resultn + row * M + col) > threshold){*(resultn + row * M + col) = 255;}
      else {*(resultn + row * M + col) = 0;}                
       }
      }

 /* each process sends its result to rank 0 process */
 MPI_Send(resultn, (M/p)*M, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);

  /* each process frees memory */
delete [] buffn;
delete[]buffn1;
delete [] resultn;
delete[] g_x;
delete[] g_y;

 }

 MPI_Finalize();

return 0;

}
  
