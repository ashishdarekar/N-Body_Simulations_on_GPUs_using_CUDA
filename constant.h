#ifndef CONSTANT_H_
#define CONSTANT_H_

#include <cuda.h>
#include <cuda_runtime.h>

#define S 1e-9f //Softening constant 
#define G 6.67408e-11 // The gravitational constant
// #define G 5 // The gravitational constant
// int numBodies = (1024 * 32);
//#define numBodies (1024 * 1024)
const float dt = (3*32*1024);   // time step for disc AU
//const float dt = 50;      //For normal simulation 
const int numIters = 20000;  // max simulation iterations
// const int cuda_blocksize = 128;  // Need to change and compare results
float totalTime = 0.0;
//float EPS2 = 1e-6f;
// #define cuda_blocksize 128
#define SOLAR_MASS 2.0e30  // in kg
#define EXTRA_MASS 1.5 // 0.02 Disk mask as a portion of center star/black hole mass
#define SYSTEM_SIZE 3.5    // Farthest particles in AU
#define SYSTEM_THICKNESS 0.08  //  Thickness in AU
#define INNER_BOUND 0.5    // Closest particles to center in AU
#define WIDTH	1024 // Image render width
#define HEIGHT	1024 // Image render height
#define BINARY_SEPARATION 0.07 // AU (only applies when binary code uncommented)
#define PI      3.14159265358979323846 
#define TO_METERS 1.496e11 // Meters in an AU
#define SOFTENING (0.015*TO_METERS) // Softens particles interactions at close distances
#define SOFTENINGSQ (SOFTENING * SOFTENING)
#define TIME_STEP (3*32*1024) //(1*128*1024) Simulated time between integration steps, in seconds
#define VTKWRITE 1
#define LOOP_UNROLL 4
#define EPS2 0.025
#define TWO_STAR 1
#define TWO_STAR_ORBIT 0.2


//Body Structure containing 7 variables
typedef struct 
{ 
  float m, x, y, z, vx, vy, vz; 
} Body;

//Body Structure containing cuda variables float4
typedef struct 
{ 
  float4 *position; //x,y,z positions w is for mass
  float4 *velocity; //vx,vy,vz velocities
} cudaBody;


#endif /* CONSTANTS_H_ */
