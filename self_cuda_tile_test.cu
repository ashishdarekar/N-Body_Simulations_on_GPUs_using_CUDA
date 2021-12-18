// implementation of CUDA code for N Body Simulations written by Ashish and Shubham 
#include <math.h>
#include <ostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda.h>
#include <random>
#include <fenv.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <sys/stat.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <string>
#include "constant.h"


// cuda error checking
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
    if (code != cudaSuccess) 
    {
        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line); 
        if (abort) exit(code);
    }
}

// Function to calculate the force and hence acceleration [24 FLOPS]
__device__ float3 acceleration(float4 posi, float4 posj, float3 F){
            
            // r_ij  [6 FLOPS]
            float dx = (posj.x - posi.x)*TO_METERS;
            float dy = (posj.y - posi.y)*TO_METERS;
            float dz = (posj.z - posi.z)*TO_METERS;
            
            // [6 Flops]
            float distSqr = dx*dx + dy*dy + dz*dz + (SOFTENINGSQ); //Numerator of the equation
            
            // invDistCube =1/distSqr^(3/2)  [4 FLOPS (2 mul, 1 sqrt, 1 inv)]
            float invDist = 1.0f / sqrtf(distSqr);
            float invDist3 = invDist * invDist * invDist; 

            //Sumation of the forces on i particle by j particle (check addition of masses is fine)
            //[ 2 + 6 = 8FLOPS]
            float gdist = G * (posj.w) * invDist3 ;

            F.x += gdist * dx ;
            F.y += gdist * dy ;
            F.z += gdist * dz ;
            return F;
}

//Kernel for Force calculation on each body. 
//FOR DISC in AU
//Tile computation 
__global__ void discbodyinteraction(float4 *pos, float4 *vel, float dt, int num_body, int cuda_Blocksize)
{
  //indexing in CUDA to capture correct body
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  
  extern __shared__ float4 shposition[]; //making shared memory of that threadblock
  if(i < num_body)
  {
    
    for (int tile = 0; tile < gridDim.x; tile++) //calculations of tiles (in a thradbloack) in a sequestial way
    { 
        //Tile preparation and copy of a data from global memory to share memory of threadblock
        float3 F;
        F.x = 0.0f;
        F.y = 0.0f;
        F.z = 0.0f;
        // float Fx = 0.0f; 
        // float Fy = 0.0f; 
        // float Fz = 0.0f;
        shposition[threadIdx.x] = pos[tile * blockDim.x + threadIdx.x];
        __syncthreads(); //To make sure all threads in a tile are copied data propely from Global memory before starting calculation
          
        //Same force calculation as that in serial
        for (int j = 0; j < cuda_Blocksize; ) 
        {
          F = acceleration(pos[i], shposition[j], F);
          j += 1;

          #if LOOP_UNROLL > 1
          F = acceleration(pos[i], shposition[j], F);
          j += 1;
          // printf("\nLoop unroll for 1!\n");
          #endif

          #if LOOP_UNROLL > 2
          F = acceleration(pos[i], shposition[j], F);
          j += 1;
          F = acceleration(pos[i], shposition[j], F);
          j += 1;
          // printf("\nLoop unroll for 2!\n");
          #endif

          #if LOOP_UNROLL > 4
          F = acceleration(pos[i], shposition[j], F);
          j += 1;
          F = acceleration(pos[i], shposition[j], F);
          j += 1;
          F = acceleration(pos[i], shposition[j], F);
          j += 1;
          F = acceleration(pos[i], shposition[j], F);
          j += 1;
          // printf("\nLoop unroll for 4!\n");
          #endif
        }
        __syncthreads(); //To make sure all threads in a tile are done with the calculation before proceed with the velocity updation

        //Acceleration and new velocitites of particle i [6 FLOPS]
        vel[i].x += ( (dt*F.x) ); 
        vel[i].y += ( (dt*F.y) ); 
        vel[i].z += ( (dt*F.z));
    }
    //Position update [1+6 = 7 FLOPS]
      pos[i].x += ( (vel[i].x)*dt / TO_METERS );
      pos[i].y += ( (vel[i].y)*dt / TO_METERS );
      pos[i].z += ( (vel[i].z)*dt / TO_METERS );
  }
}

__global__ void initializer(float4 *pos, float4 *vel, const int num_body){
  curandState state;
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  #if TWO_STAR
  if (i==0){
    pos[i].x = TWO_STAR_ORBIT;
    pos[i].y = 0.0;
    pos[i].z = 0.0;
    pos[i].w = 100*SOLAR_MASS;
    vel[i].x = 0.0;
    vel[i].y = -pow(((G*(SOLAR_MASS + EXTRA_MASS*SOLAR_MASS))/ (TWO_STAR_ORBIT*TO_METERS)), 0.5);
    vel[i].z = 0.0;
  }
  else if (i==1){
    pos[i].x = -TWO_STAR_ORBIT;
    pos[i].y = 0.0;
    pos[i].z = 0.0;
    pos[i].w = 100*SOLAR_MASS;
    vel[i].x = 0.0;
    vel[i].y = pow(((G*(SOLAR_MASS + EXTRA_MASS*SOLAR_MASS))/ (TWO_STAR_ORBIT*TO_METERS)), 0.5);
    vel[i].z = 0.0;
  }
  #else
  if (i==0){
    pos[i].x = 0.0;
    pos[i].y = 0.0;
    pos[i].z = 0.0;
    pos[i].w = SOLAR_MASS;
    vel[i].x = 0.0;
    vel[i].y = 0.0;
    vel[i].z = 0.0;
  }
  #endif
  else{
    float factor;
    
    #if TWO_STAR
    factor = 1.0;
    #else
    factor = 1.0;
    #endif
    
    float angle;
    float radius, randRadius;
    float velocity, randHeight;
    velocity = 0.67*sqrt((G*SOLAR_MASS)/(4*BINARY_SEPARATION * TO_METERS));
    // reference article https://stackoverflow.com/questions/18501081/generating-random-number-within-cuda-kernel-in-a-varying-range
    curand_init(clock64(), i, 0, &state);
    angle = curand_uniform(&state)*(2.0*PI - 0.0);
    curand_init(clock64(), i, 0, &state);
    randRadius = curand_uniform(&state)*(SYSTEM_SIZE - INNER_BOUND) + INNER_BOUND;
    curand_init(clock64(), i, 0, &state);
    randHeight = curand_uniform(&state)*(SYSTEM_THICKNESS);

    radius = sqrtf(SYSTEM_SIZE)*sqrtf(randRadius);
		velocity = pow((G*((200*SOLAR_MASS)+1*(EXTRA_MASS*SOLAR_MASS))/ (radius*TO_METERS)), 0.5);

    pos[i].x = radius*cos(angle);
    pos[i].y = radius*sin(angle);
    pos[i].z = randHeight-SYSTEM_THICKNESS/2.0;
    pos[i].w = (EXTRA_MASS*SOLAR_MASS)/num_body;
    vel[i].x = velocity*sin(angle);
    vel[i].y = -velocity*cos(angle);
    vel[i].z = 0.0;
  }
}

//Visualisation using VTk
void write_vtkFile(const char *szProblem, int timeStepNumber, int num_body, cudaBody *p ) 
{
  char szFileName[80];
  FILE *fp=NULL;
  sprintf( szFileName, "%s.%i.vtk", szProblem, timeStepNumber );
  fp = fopen( szFileName, "w");
  if( fp == NULL )		       
  {
    char szBuff[80];
    sprintf( szBuff, "Failed to open %s", szFileName );
    //ERROR( szBuff );
    return;
  }

  // Write VTK Header
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"generated for CUDA Seminar output (written by Ashish and Shubham) \n");
  fprintf(fp,"ASCII\n");
  fprintf(fp,"DATASET POLYDATA\n");
  fprintf(fp,"POINTS %i float\n", num_body);

  // Write VTK point Coordinator
  for( int i = 0; i < num_body; i++ ) 
  {
    fprintf(fp, "%f %f %f\n",(p->position[i].x), (p->position[i].y), (p->position[i].z));
  }
	
  if( fclose(fp) )
  {
    char szBuff[80];
    sprintf( szBuff, "Failed to close %s", szFileName );
    //ERROR( szBuff );
  }
}

//Function to calculate the Performance in GFLOPS
double computePerf(int cuda_numblocks, float timeseconds, const int iterations, int num_body)
{
    int flopsPerInteraction=0;
    flopsPerInteraction = ( num_body * num_body * 24 ) + ( num_body *  (num_body/cuda_numblocks) * 6 ) + ( num_body * 7);
    
    double total_operations = iterations * (double)flopsPerInteraction ; 
    double gflops = 1e-9 * ((double)total_operations) / timeseconds  ;
    
    return gflops;
}


int main(int argc, char* argv[])
{
  int numBodies;
  int cuda_blocksize;

  if (argc>=2)
  {
    if (std::string(argv[1]) == "-s") 
    {      
      numBodies = std::stoi(argv[2]);  
    } 
    else
    {
      fprintf(stderr, "Wrong argument - %s\n", argv[1]);
      exit(EXIT_FAILURE);
    }
    if (std::string(argv[3]) == "-b"){
      cuda_blocksize = std::stoi(argv[4]);
    }
    else
    {
      fprintf(stderr, "Wrong argument - %s\n", argv[4]);
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    fprintf(stderr, "Usage:./self_cuda_tile -s <numBodies> -b <blocksize>\n");
    exit(EXIT_FAILURE);
  }
  //To save all the visualisation files in the folder
  char sol_directory[100];
  char sol_folder[100];
  struct stat st;

#if VTKWRITE

  sprintf(sol_folder,"cuda_solution_%s","cuda_seminar");
  if(stat(sol_folder,&st)==-1)
    mkdir(sol_folder,0700);

  sprintf(sol_directory,"cuda_solution_%s/sol","cuda_seminar");

#endif

  //Memory allocation on the HOST
  int size_mem = numBodies * 2 * sizeof(float4);
  float4 *h_buf = (float4*) malloc(size_mem);
  cudaBody ph = {h_buf, h_buf + numBodies }; //2 arrays, positions and velocity arrays, so ph pointing to the h_buf

  //Memory allocation on the DEVICE
  float4 *d_buf;
  cudaMalloc((void**) &d_buf, size_mem );
  cudaBody pd = { d_buf, d_buf + numBodies }; 

  int cuda_numblocks = (numBodies + cuda_blocksize - 1) / cuda_blocksize;  

  //Time Iterations 
  clock_t startTime = clock();
  dim3 grid = dim3(cuda_numblocks, 1, 1);
  dim3 threads = dim3(cuda_blocksize, 1, 1);

  initializer<<<grid, threads>>>(pd.position, pd.velocity, numBodies);
  
  cudaMemcpy(h_buf, d_buf, size_mem , cudaMemcpyDeviceToHost);
  write_vtkFile(sol_directory, 0, numBodies, &ph);

  // for timining
  float time;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);//Time Begin
  
  for (int iter=1; iter<=numIters; iter++)
  {
    //Copy data from Host to Device, buffer copying 
    discbodyinteraction<<<grid, threads, cuda_blocksize*sizeof(float4)>>>( pd.position, pd.velocity, dt, numBodies, cuda_blocksize);
  } 

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);
  printf("\nElapsed time for Computation:  %10.5f sec \n", time/1000.0f);
  double gflops=0.0;
  gflops = computePerf(cuda_numblocks, time/1000.0f, numIters, numBodies);
  
  printf("\nCUDA N-Body (%d bodies): %f GFLOP/s\n", numBodies, gflops);

  cudaMemcpy(h_buf, d_buf, size_mem , cudaMemcpyDeviceToHost);
  write_vtkFile(sol_directory, numIters, numBodies, &ph);

  free(h_buf);
  cudaFree(d_buf);
}

