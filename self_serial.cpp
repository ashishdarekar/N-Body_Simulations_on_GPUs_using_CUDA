// Serial Version of Particle Partcile interaction written by "Ashish" and "Shubham"

#include <math.h>
#include <ostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <random>
#include <fenv.h>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string>
#include "constant.h"
#include <omp.h>
//#include <cuda.h>
//#include <cuda_runtime.h>
//#include <curand.h>
//#include <curand_kernel.h>


//FOR DISC in AU
//function to initialize all 7 parameters mass, positions and velocities
void discinitilize(Body *p, const int num_body)
{
  using std::uniform_real_distribution;
	uniform_real_distribution<double> randAngle (0.0, 200.0*PI);
	uniform_real_distribution<double> randRadius (INNER_BOUND, SYSTEM_SIZE);
	uniform_real_distribution<double> randHeight (0.0, SYSTEM_THICKNESS);
	std::default_random_engine gen (0);
	double angle;
	double radius;
	double velocity;
	Body *current;

	//STARS
	velocity = 0.67*sqrt((G*SOLAR_MASS)/(4*BINARY_SEPARATION * TO_METERS));
	//STAR 1 is SUN
	current = &p[0];
	current->x = 0.0;///-BINARY_SEPARATION;
	current->y = 0.0;
	current->z = 0.0;
	current->vx = 0.0;
	current->vy = 0.0;//velocity;
	current->vz = 0.0;
	current->m = SOLAR_MASS;


	///STARTS AT NUMBER OF STARS///
	double totalExtraMass = 0.0;
	#pragma omp parallel for schedule(static, 8) private(angle, radius, velocity) reduction(+: totalExtraMass)
	for (int index=1; index<num_body; index++)
	{
		angle = randAngle(gen);
		radius = sqrt(SYSTEM_SIZE)*sqrt(randRadius(gen));
		velocity = pow(((G*(SOLAR_MASS+((radius-INNER_BOUND)/SYSTEM_SIZE)*EXTRA_MASS*SOLAR_MASS))
					  	  	  	  	  / (radius*TO_METERS)), 0.5);
		current = &p[index];
		current->x =  radius*cos(angle);
		current->y =  radius*sin(angle);
		current->z =  randHeight(gen)-SYSTEM_THICKNESS/2;
		current->vx =  velocity*sin(angle);
		current->vy = -velocity*cos(angle);
		current->vz =  0.0;
		current->m = (EXTRA_MASS*SOLAR_MASS)/num_body;
		totalExtraMass += (EXTRA_MASS*SOLAR_MASS)/num_body;
	}
	// std::cout << "\nTotal Disk Mass: " << totalExtraMass;
	// std::cout << "\nEach Particle weight: " << (EXTRA_MASS*SOLAR_MASS)/num_body
	// 		  << "\n______________________________\n";
}


//FOR DISC in AU
//Force calculation for Disc bodies
void discbodyinteraction(Body *p, float dt, const int num_body)
{
  
  for (int i = 0; i < num_body; i++) 
  { 
    float Fx = 0.0f; 
    float Fy = 0.0f; 
    float Fz = 0.0f;
    float dx, dy, dz, invDist, invDist3;
	  double distSqr;	

    #pragma omp parallel for simd schedule(static, 8) private(dx, dy, dz, invDist, invDist3, distSqr) reduction(+: Fx, Fy, Fz)
    for (int j = 0; j < num_body; j++) 
    {
      dx = (p[j].x - p[i].x)*TO_METERS;
      dy = (p[j].y - p[i].y)*TO_METERS;
      dz = (p[j].z - p[i].z)*TO_METERS;
      distSqr = dx*dx + dy*dy + dz*dz + (SOFTENING * SOFTENING); //Numerator of the equation
      invDist = 1.0f / sqrt(distSqr);
      invDist3 = invDist * invDist * invDist; 

      //Sumation of the forces on i particle by j particle (check addition of masses is fine)
      Fx += G * (p[i].m) * (p[j].m) * dx * invDist3;
      Fy += G * (p[i].m) * (p[j].m) * dy * invDist3;
      Fz += G * (p[i].m) * (p[j].m) * dz * invDist3;
    }

    //Acceleration and new velocitites of particle i
    p[i].vx += ( (dt* Fx) / (p[i].m) ); 
    p[i].vy += ( (dt*Fy) / (p[i].m) ); 
    p[i].vz += ( (dt*Fz) / (p[i].m) );
    }

}

//FOR DISC in AU
//function for updating the postion according to the Euler time stepping
void discupdateposition(Body *p, float dt, const int num_body)
{
  #pragma omp parallel for simd schedule(static, 8)
  for (int i = 0 ; i < num_body; i++) 
  { 
      p[i].x += p[i].vx*dt/TO_METERS;
      p[i].y += p[i].vy*dt/TO_METERS;
      p[i].z += p[i].vz*dt/TO_METERS;
  }
}

//Visualisation using VTk
void write_vtkFile(const char *szProblem, int timeStepNumber, int num_body, Body *p ) 
{
  int i,j;
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
  fprintf(fp,"generated for CUDA Seminar output (written by Ashish) \n");
  fprintf(fp,"ASCII\n");
  fprintf(fp,"DATASET POLYDATA\n");
  fprintf(fp,"POINTS %i float\n", num_body);

  // Write VTK point Coordinator
  for( int i = 0; i < num_body; i++ ) 
  {
    fprintf(fp, "%f %f %f\n",(p[i].x), (p[i].y), (p[i].z) );
  }
	
  if( fclose(fp) )
  {
    char szBuff[80];
    sprintf( szBuff, "Failed to close %s", szFileName );
    //ERROR( szBuff );
  }
}


int main(int argc, char* argv[])
{
  int numBodies;

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
  }
  else
  {
    fprintf(stderr, "Usage:./self_serial -s <numBodies>\n");
    exit(EXIT_FAILURE);
  }

  //To save all the visualisation files in the folder
  char sol_directory[100];
  char sol_folder[100];
  struct stat st;
  sprintf(sol_folder,"serial_solution_%s","cuda_seminar");
  if(stat(sol_folder,&st)==-1)
    mkdir(sol_folder,0700);

  sprintf(sol_directory,"serial_solution_%s/sol","cuda_seminar");

  //Array of all bodies
  Body *p = (Body*) malloc(numBodies * sizeof(Body));

  //Initialisation of all parameters realted to body
  discinitilize(p, numBodies);

  //Initial Data
  //write_vtkFile(sol_directory,0,numBodies,p);

  //Time Iterations 
#ifdef _OPENMP
  double before;
  if(omp_get_thread_num()==0){
    before = omp_get_wtime();
  }
#else
  clock_t before= clock();
#endif

  int msec=0;

  for (int iter=1; iter<=numIters; iter++)
  {
    //Calculate all pairwise interactions
    discbodyinteraction(p, dt, numBodies);

    //Update positons
    discupdateposition(p, dt, numBodies);

    /*
    //VTK writing
    if(iter % 100 == 0 ) 
    {
      write_vtkFile(sol_directory,iter,numBodies,p);
      printf("%d iterations done. Printing files. \n",iter);
    }
    //printf("Iteration %d done\n", iter);
    */
  }

  //End time Iterations
  //write_vtkFile(sol_directory,numIters,numBodies,p);
#if _OPENMP
  double difference;
  if(omp_get_thread_num()==0){
    difference = omp_get_wtime() - before;
    printf("Total time for %d Iterations : %f seconds\n", numIters, difference);
  }
#else
  clock_t difference = clock() -before;
  msec = difference * 1000 / CLOCKS_PER_SEC;
  printf("Total time for %d Iterations : %d seconds\n", numIters, msec/1000);
#endif
    
  free(p);
   
}
