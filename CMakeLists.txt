cmake_minimum_required(VERSION 3.0.0)

project(SEWAS)

# version number
set(SEWAS_VERSION_MAJOR 1)
set(SEWAS_VERSION_MINOR 0)

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules" CACHE STRING "")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are None, Debug, Release, RelWithDebInfo and MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

add_compile_options(-mcx16)
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  add_compile_options(-march=native -mtune=native)
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-funroll-loops>)
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-march=native>)
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-mtune=native>)
  # add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DEIGEN_ENABLE_AVX512>)
endif()

option(COLLECT_STATS "Collect statistics during execution (default ON)" ON)
if (COLLECT_STATS)
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DCOLLECT_STATS>)

  option(SHOW_PER_THREAD_STATS "Enable showing statistics at thread-level (default OFF)" OFF)
  if(SHOW_PER_THREAD_STATS)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DSHOW_PER_THREAD_STATS>)
  endif(SHOW_PER_THREAD_STATS)
endif(COLLECT_STATS)

option(PLOT_DATA "Enable ploting the data of a tile object (default OFF)" OFF)
if (PLOT_DATA)
  # Adjust compile options accordingly
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DUSE_MATPLOTLIB>)
endif(PLOT_DATA)

option(VISUALIZE_EXECUTION "Visualize the tasks execution (default OFF)" OFF)
if(VISUALIZE_EXECUTION)
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DUSE_VTK>)
  option(ENABLE_VELOCITY_RENDERING "Visualize the execution of velocity tasks (default OFF)" ON)
  option(ENABLE_STRESS_RENDERING "Visualize the execution of stress tasks (default OFF)" OFF)
  option(ENABLE_CLUSTER_RENDERING "Visualize the cluster activity during tasks execution (default OFF)" OFF)

  # Adjust compile options accordingly
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DOPENGL2_SUPPORTED>) ## Not working for Raspberry PI 3
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DVISUALIZE_EXECUTION>)
  if (ENABLE_VELOCITY_RENDERING)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DENABLE_VELOCITY_RENDERING>)
  endif(ENABLE_VELOCITY_RENDERING)
  if (ENABLE_STRESS_RENDERING)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DENABLE_STRESS_RENDERING>)
  endif(ENABLE_STRESS_RENDERING)
  if (ENABLE_CLUSTER_RENDERING)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DENABLE_CLUSTER_RENDERING>)
  endif(ENABLE_CLUSTER_RENDERING)
endif(VISUALIZE_EXECUTION)

option(USE_BLOCKWISE_FDO "Use block operations for evaluating the finite difference operator (default ON)" ON)
if (USE_BLOCKWISE_FDO)
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DBLOCKWISE_FDO>)
endif(USE_BLOCKWISE_FDO)

option(VECTORIZE_COMPUTATIONS "Explicit vectorization of the computation (default ON)" ON)
if(VECTORIZE_COMPUTATIONS)
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DVECTORIZE_COMPUTATIONS>)

  option(EIGEN_VECTORIZATION "Use Eigen as the vectorization engine (default ON)" ON)
  if (EIGEN_VECTORIZATION)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DEIGEN_VECTORIZATION>)
  endif(EIGEN_VECTORIZATION)

  option(BOOST_SIMD_VECTORIZATION "Use Boost.SIMD as the vectorization engine (default OFF)" OFF)
  if (BOOST_SIMD_VECTORIZATION)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DBOOST_SIMD_VECTORIZATION>)
  endif(BOOST_SIMD_VECTORIZATION)
endif(VECTORIZE_COMPUTATIONS)

# ----------------------------------------------
# Build configuration completed
# ----------------------------------------------

option(SEWAS_WITH_PARSEC "Use PaRSEC for task-scheduling (default ON)" ON)
if(SEWAS_WITH_PARSEC)
  add_compile_options(-DSEWAS_WITH_PARSEC)

  # Compilation of the JDF
  include(RulesJDF)
  find_program(parsec_ptgpp_EXE parsec_ptgpp)

  set(PARSEC_DEP_MGT "HASH_TABLE" CACHE STRING "PaRSEC dependency management (INDEX_ARRAY|HASH_TABLE)")
  mark_as_advanced(FORCE PARSEC_DEP_MGT)
  if ("${PARSEC_DEP_MGT}" STREQUAL "INDEX_ARRAY")
    set(PARSEC_PTGPP_CFLAGS "-M;index-array") # Workaround suggested by George to avoid sporadic non-termination
  endif()

  jdf_rules(JDF2C_SOURCES ${CMAKE_SOURCE_DIR}/src/runtimes/parsec/sewas.jdf)
  include_directories(${CMAKE_CURRENT_BINARY_DIR})
endif(SEWAS_WITH_PARSEC)

set(SOURCES
  src/main.cxx
  src/CartesianMesh3D.cxx
  src/DataSet.cxx
  src/ExternalSource.cxx
  src/HaloManager.cxx
  src/LinearSeismicWaveModel.cxx
  src/SEWASParameterManager.cxx
  src/SEWASSequential.cxx
  src/SEWASPaRSEC.cxx
  src/Mesh3DPartitioning.cxx
  src/VisualizationManager.cxx
  src/MetricsManager.cxx)

include_directories(include)

add_executable(sewas ${SOURCES} ${JDF2C_SOURCES})
set_property(TARGET sewas PROPERTY CXX_STANDARD 14)

if(SEWAS_WITH_PARSEC)
  # PaRSEC
  find_package(PARSEC 2.0 REQUIRED)
  include_directories("${PARSEC_INCLUDE_DIRS}")
  target_link_libraries(sewas ${PARSEC_LIBRARIES})
endif(SEWAS_WITH_PARSEC)

# Boost
find_package(Boost 1.62 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
target_link_libraries(sewas ${Boost_LIBRARIES})

# Eigen
find_package(Eigen3 3.3.3 REQUIRED)
include_directories("${EIGEN3_INCLUDE_DIR}")

# Boost.SIMD
if (BOOST_SIMD_VECTORIZATION)
  find_package(Boost.SIMD 4.17 REQUIRED)
  include_directories("${Boost.SIMD_INCLUDE_DIRS}")
endif(BOOST_SIMD_VECTORIZATION)

# We rely on Intel TBB concurrent containers for collecting metrics
if (COLLECT_STATS)
  find_package(TBB 4.3 REQUIRED)
  include_directories("${TBB_INCLUDE_DIRS}")
  target_link_libraries(sewas ${TBB_LIBRARIES})
endif(COLLECT_STATS)

# The package matplotlib-cpp is used to plot data -> https://github.com/lava/matplotlib-cpp
# TODO Add an alternative for using xplot -> https://github.com/QuantStack/xplot
if (PLOT_DATA)
  include_directories("${CMAKE_SOURCE_DIR}/thirdparty/matplotlibcpp")
  find_package(PythonLibs 2.7)
  target_include_directories(sewas PRIVATE ${PYTHON_INCLUDE_DIRS})
  target_link_libraries(sewas ${PYTHON_LIBRARIES})
endif(PLOT_DATA)

# VTK is used to visualize the tasks execution
if (VISUALIZE_EXECUTION)
  find_package(VTK 8.1 REQUIRED)
  include_directories("${VTK_INCLUDE_DIRS}")
  target_link_libraries(sewas ${VTK_LIBRARIES})
endif(VISUALIZE_EXECUTION)
