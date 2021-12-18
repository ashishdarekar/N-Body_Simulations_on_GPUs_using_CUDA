**Fast N-Body Simulation with CUDA**

Seminar Work : "Parallelization of Physics calculations on GPUs with CUDA"
By Ashish Darekar and Shubham Khatri
Technische Universität München, Fakultät für Physik

To use the code, you will need The NVIDIA CUDA Toolkit version 
1.0 or later, available at http://developer.nvidia.com/cuda

**Building the Code**   
Linux:
1. Go to the cuda_seminar_nbody_2020 directory 
2. Go to "constant.h" : To change the problem input parameters like Number of iterations, VTKWRITE flag (Default parameters are provided)
3. make clean
4. make 

**Running the code for Visualization**
(With Visualization of Results)
(VTKWRITE 1)

There are three execution codes: (Command options must provide)
1. CPU Serial Code:           ./serial -s 1024
2. CPU Code with OpenMP:      ./serial_omp -s 1024
3. GPU Code:                  ./self_cuda_tile -s <Number of bodies> -b <Block size>
        eg.:                  ./self_cuda_tile -s 16384 -b 128

Results: VTK files stored in "serial_solution_cuda_seminar" for Serial Code and Serial With OpenMP
                             "cuda_solution_cuda_seminar" for CUDA code
         Can be locally visualize ParaView
      
**Running this code for Performance Benchmarking**
Performance Benchmarking: Python scripts are provided (exec.sh, serial_exec.sh, serial_omp.sh)
(VTKWRITE 0)

Performance Benchmarking with NVIDIA code: Python script is provided (perf.sh)
(VTKWRITE 0)

         
 
