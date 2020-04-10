![C/C++ CI](https://github.com/aneoconsulting/SeWaS/workflows/C/C++%20CI/badge.svg?branch=master)

# General Description of SeWaS
This document gives an overview of [SeWaS](https://github.com/aneoconsulting/SeWaS), a modern production-class Seismic Wave Simulator. This simulator is built from scratch using C++14 and makes extensive use of state-of-the-art libraries such as [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) and [Boost](http://www.boost.org/). In addition, it embeds advanced techniques for developing HPC applications on top of emerging architectures, such as cache-blocking and vectorization. In the following, we give a macroscopic description of the underlying algorithm, data structures and the parallelization scheme.

# User Interface
+ Currently, the application requires a cubic computational domain. The domain description and mechanical parameters are defined in a JSON file as following:

```json
{
    "Model" : "Dummy",
    "tmax" : 1.6,
    "dt" : 0.008,
    "lx" : 60000,
    "ly" : 60000,
    "lz" : 20500,
    "ds" : 100,
    "Layers" :
    [
        {
            "Name" : "Basin 1",
            "Size" :
            {
                "X" :
                {
                    "start" : 0,
                    "end" : 60000
                },
                "Y" :
                {
                    "start" : 0,
                    "end" : 60000
                },
                "Z" :
                {
                    "start" : 0,
                    "end" : 300
                }
            },
            "Vp" : 600,
            "Vs" : 300,
            "rho" : 1700,
            "Q" : 70
        },
        {
            "Name" : "Basin 2",
            "Size" :
            {
                "X" :
                {
                    "start" : 0,
                    "end" : 60000
                },
                "Y" :
                {
                    "start" : 0,
                    "end" : 60000
                },
                "Z" :
                {
                    "start" : 300,
                    "end" : 750
                }
            },
            "Vp" : 1870,
            "Vs" : 1000,
            "rho" : 2000,
            "Q" : 100
        },
        {
            "Name" : "Basin 3",
            "Size" :
            {
                "X" :
                {
                    "start" : 0,
                    "end" : 60000
                },
                "Y" :
                {
                    "start" : 0,
                    "end" : 60000
                },
                "Z" :
                {
                    "start" : 750,
                    "end" : 1500
                }
            },
            "Vp" : 2800,
            "Vs" : 1500,
            "rho" : 2300,
            "Q" : 100
        },
        {
            "Name" : "Bedrock 1",
            "Size" :
            {
                "X" :
                {
                    "start" : 0,
                    "end" : 60000
                },
                "Y" :
                {
                    "start" : 0,
                    "end" : 60000
                },
                "Z" :
                {
                    "start" : 1500,
                    "end" : 5000
                }
            },
            "Vp" : 5000,
            "Vs" : 2890,
            "rho" : 2700,
            "Q" : 200
        },
        {
            "Name" : "Bedrock 2",
            "Size" :
            {
                "X" :
                {
                    "start" : 0,
                    "end" : 60000
                },
                "Y" :
                {
                    "start" : 0,
                    "end" : 60000
                },
                "Z" :
                {
                    "start" : 5000,
                    "end" : 20000
                }
            },
            "Vp" : 6000,
            "Vs" : 3465,
            "rho" : 2700,
            "Q" : 250
        }
    ],
    "Sources" :
    [
        {
            "Name" : "A",
            "xs" : 1000,
            "ys" : 1000,
            "zs" : 1000
        }
    ]
}
```
+ The application is launched from command line as following:

```sh
    $ ./sewas --cx 40 --cy 40 --cz 40 --P 1 --Q 1 --R 1 --nthreads 3 --dfile=../data/input/dummy.json
```

where the arguments are detailed below:
```sh
Seismic Wave Simulator (SeWaS)

Allowed options:
  -h [ --help ]         Produce help message
  -x [ --cx ] arg       Block size along x-axis
  -y [ --cy ] arg       Block size along y-axis
  -z [ --cz ] arg       Block size along z-axis
  -P [ --P ] arg        Number of MPI processes along x-axis
  -Q [ --Q ] arg        Number of MPI processes along y-axis
  -R [ --R ] arg        Number of MPI processes along z-axis
  -T [ --nthreads ] arg Number of threads per MPI process
  -D [ --dfile ] arg    Path to the Json file containing mechanical properties
                        and kinematic source parameter
  -C [ --config ] arg   Configuration file
```

+ An alternative way to launch the application consists to use a configuration file containing all required input parameters. The previous command line can then be simplified to:
``` sh
$ ./sewas -C dummy.cfg
```

where dummy.cfg contains:
```sh
cx=40
cy=40
cz=40
P=1
Q=1
R=1
nthreads=3
dfile=../data/input/dummy.json
```

# Algorithm
The application implements the numerical solution of the linear elastodynamic equation. The discretization follows a staggered 4th order finite difference scheme in space and 2nd order in time. The numerical schemes considered are detailed in [Dupros2010](https://tel.archives-ouvertes.fr/tel-00580411/).

# Data Structures
The application is developed using Eigen's containers. The spatial domain is decomposed into a collection of contiguous 3D tiles of cells. Each tile object is an instance of the SpatialBlockField class, a multi-dimensional
container mapped on a 1D `Eigen::Array` storing cx*cy*cz cells. In the following, we show the descriptions of the Velocity, Stress and SourceForce fields.

ii,jj,kk variables are used to represent a tile coordinates, whereas i,j,k represent a cell coordinates. The time-step is denoted ts.

## Velocity
```c++
  // (x|y|z)(ii, jj, kk)(i, j, k)
  typedef Eigen::Array<SpatialField, DIM, 1, Eigen::ColMajor> Velocity;
```

Let us consider an object `V` of type Velocity. `Vx=V(x)(ii,jj,kk)` is a 3D spatial block of cx\*cy\*cz cells. Such an approach allows to perform *block-wise operations*:

```c++
  Vx=42.; // All elements of Vx will be initialized to 42.
  Vx+=U;  // U is another container of the same type as Vx
```

Alternatively, one can go through all elements of `Vx` as following:

```c++
  for (int i=0; i<cx; i++)
    for (int j=0; j<cy; j++)
      for (int k=0; k<cz; k++)
        Vx(i,j,k)=42.;
```

With this *block-wise* formulation, the x-component of the velocity field is computed as following:

```c++
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        vX(i,j)+=(fdo_.apply<L2S::X>(sigmaXX, i, j, L2S::DK)
                + fdo_.apply<L2S::Y>(sigmaXY, i, j, L2S::DK)
                + fdo_.apply<L2S::Z>(sigmaXZ, i, j, L2S::DK))*dt*bx(i,j);
      }
    }
```

where `fdo_` is an instance of the 4th order *finite difference operator* acting on a 1D spatial block of cells (for instance `vX(i,j)` is a 1D vector along z-axis).

An important point to note here is that the right-hand side, of the above mentioned velocity-update instruction, is evaluated thanks to the C++ idiom *Expression Templates*. Therefore, only one new `SpatialBlockField1D` object is created for the whole expression.

## Stress
```c++
  // (xx|yy|zz|xy|xz|yz)(ii, jj, kk)(i, j, k)
  typedef Eigen::Array<SpatialField, NB_STRESS_FIELD_COMPONENTS, 1, Eigen::ColMajor> StressField;
```

Similarly, `Sxx=S(xx)(ii,jj,kk)` is a 3D spatial block of cx\*cy\*cz cells, and all previously listed operations are applicable.

## Velocity source
```c++
// (is, js, ks) are the coordinates of a cell-source
v_(L2S::X)(ii,jj,kk)(is,js,ks)-=a*((ts-2)*dt-1.2/7.)*exp(-M_PI*M_PI*7.*7.*((ts-2)*dt-1.2/7.)*((ts-2)*dt-1.2/7.))*dt/rho(ii,jj,kk)(is,js,ks);
```

# Parallelization Scheme
The algorithm is implemented using a *task-based approach* and the parallelization is handled by the [PaRSEC](http://icl.utk.edu/parsec/) runtime system. We defined 6 computational tasks within the dataflow:

```c++
ComputeVelocity
ExtractVelocityHalo
UpdateVelocity
ComputeStress
ExtractStressHalo
UpdateStress
```


Two visualization tasks are added to enable real-time monitoring of the ongoing simulation:

```c++
DisplayVelocity
DisplayStress
```

Those visualization tasks interact with [VTK](http://www.vtk.org/) actors to render the scene.

# Building the application
The application has been tested on Debian 9 and CentOS. The following steps are suitable for a Debian-like OS. Our experiments were carried out using Intel compiler 18, but according to our recent experiments using gcc 7.3 there is no significant performance differences between these two compilers. So all the following experiments can be obtained using gcc.

## Dependencies
+ Boost 1.62 (program-options module)
```sh
$ sudo apt-get install libboost-program-options-dev
```

+ OpenMPI 2.0.2
```sh
$ sudo apt-get install libopenmpi-dev openmpi-common
```

+ PaRSEC 2.0
```sh
$ git clone https://bitbucket.org/icldistcomp/parsec.git
$ cd parsec
$ cmake .. -DPARSEC_WITH_DEVEL_HEADERS=1 -DPARSEC_GPU_WITH_CUDA=0
$ make install
```

+ Eigen 3.3.3
```sh
$ sudo apt-get install libeigen3-dev
```

+ Intel TBB 4.3
```sh
$ sudo apt-get install libtbb-dev
```
     
## SeWaS build
```sh
$ export SEWAS_ROOT=/path/to/SeismicWaveSimulator
$ cd $SEWAS_ROOT/build
$ cmake ../src -DCMAKE_BUILD_TYPE=Release
$ make
```

At this stage the `build/` directory should contain an executable `sewas`

/!\ The results presented in our paper are obtained on a cluster of Intel KNL processors. Hence it is mandatory to build the SeWaS application natively on one co-processor for best performances.

# Benchmarks
Performance studies have been carried out using two test cases: TestA and TestB.
As an example the **TestA** input file is given by:
```json
{
    "Model" : "TestA",
    "tmax" : 0.8,
    "dt" : 0.008,
    "lx" : 20000,
    "ly" : 20000,
    "lz" : 10000,
    "ds" : 100,
    "Layers" :
    [
        {
            "Name" : "Basin 1",
            "Size" :
            {
                "X" :
                {
                    "start" : 0,
                    "end" : 20000
                },
                "Y" :
                {
                    "start" : 0,
                    "end" : 20000
                },
                "Z" :
                {
                    "start" : 0,
                    "end" : 10000
                }
            },
            "Vp" : 6000,
            "Vs" : 3460,
            "rho" : 2700
        }
    ],
    "Sources" :
    [
        {
            "Name" : "A",
            "xs" : 1000,
            "ys" : 1000,
            "zs" : 1000,
            "Depth" : 5000
        }
    ]
}
```

The results of the numerical computations using **TestA** are validated with the reference application Ondes3D. The statistics relative to this experiment, conducted on a dual-core laptop are given below.

```sh
salli@jigawal:~/Projects/HPC/SeismicWaveSimulator/build$ ./sewas -C TestA.cfg 
Source : (13, 13, 13) within tile (0, 0, 0)
Seismic wave is progressing...
I am 0 and I completed the simulation...
++++++++++++++++++++++++++++++++++++++++++
 Statistics 
------------------------------------------
 SEWAS 
++++++++++++++++++++++++++++++++++++++++++
                                AVG	        MIN	        MAX	        TOT	        SD          
Global
	#Calls           : 	1
	Elapsed Time (s) : 	1.743005e+01	1.743005e+01	1.743005e+01	1.743005e+01	0.000000e+00
	CPU Time (s)     : 	3.391434e+01	3.391434e+01	3.391434e+01	3.391434e+01	0.000000e+00
Core simulation
	#Calls           : 	1
	Elapsed Time (s) : 	1.681025e+01	1.681025e+01	1.681025e+01	1.681025e+01	0.000000e+00
	CPU Time (s)     : 	3.357347e+01	3.357347e+01	3.357347e+01	3.357347e+01	0.000000e+00
Initialization
	#Calls           : 	1
	Elapsed Time (s) : 	5.559187e-01	5.559187e-01	5.559187e-01	5.559187e-01	0.000000e+00
	CPU Time (s)     : 	2.772830e-01	2.772830e-01	2.772830e-01	2.772830e-01	0.000000e+00
ComputeStress
	#Calls           : 	9504
	Elapsed Time (s) : 	2.418401e-03	8.470220e-04	1.123556e-02	2.298448e+01	2.183366e-03
	CPU Time (s)     : 	4.808716e-03	8.600000e-04	2.409800e-02	4.570203e+01	4.384567e-03
ComputeVelocity
	#Calls           : 	4752
	Elapsed Time (s) : 	1.493992e-03	1.189589e-03	3.632833e-03	7.099449e+00	5.015644e-04
	CPU Time (s)     : 	3.005864e-03	1.190000e-03	8.809000e-03	1.428387e+01	1.283688e-03
```

# Results reproducibility
For a distributed run, the P, Q and R parameters (defined in both **TestA.cfg** and **TestB.cfg** files) need to be changed according to the total number of MPI processes used, such that:
```sh
NB_MPI_PROCESSES=P*Q*R
```

We provide the best (P,Q,R) used to run our experiments on TestB (results presented in section 5 of the paper):

| NB_MPI_PROCESSES | P | Q | R |
| -----------------|---|---|---|
|                1 | 1 | 1 | 1 |
|                2 | 2 | 1 | 1 |
|                4 | 2 | 2 | 1 |
|                8 | 2 | 4 | 1 |



For example on 8 nodes, the application is launched as following:
```sh
mpirun -n 8 -pernode ./sewas --cx 40 --cy 40 --cz 100 --P 2 --Q 4 --R 1 --nthreads=64 --dfile=../data/input/TestB.json
```
