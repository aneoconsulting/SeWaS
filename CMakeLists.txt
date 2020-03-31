cmake_minimum_required(VERSION 3.10.3)

project(SEWAS)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# version number
set(SEWAS_VERSION_MAJOR 0)
set(SEWAS_VERSION_MINOR 1)
set(SEWAS_VERSION_PATCH 0)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules" CACHE STRING "")

option(SEWAS_DISTRIBUTED "Build SeWaS for distributed memory platforms (default OFF)" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are None, Debug, Release, RelWithDebInfo and MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release" OR "${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    add_compile_options(-march=native -mtune=native)
		add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-funroll-loops>)
		add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-march=native>)
		add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-mtune=native>)
	elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	endif()

	add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DEIGEN_ENABLE_AVX512>)
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  add_compile_options(-mcx16)
  add_compile_options(-m64)
endif()

option(VERBOSE "Enable collecting and displaying of traced events at runtime (default OFF)" OFF)
if (VERBOSE)
  add_compile_options(-DVERBOSE)
endif(VERBOSE)

option(BUILD_TESTING "Build the testing tree (default OFF)" OFF)
if (BUILD_TESTING)
  add_compile_options(-DBUILD_TESTING)
endif(BUILD_TESTING)

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

option(ENABLE_IO "Enable performing IO at runtime (default OFF)" OFF)
if (ENABLE_IO)
  # Adjust compile options accordingly
  add_compile_options(-DENABLE_IO)
endif(ENABLE_IO)

option(VISUALIZE_EXECUTION "Visualize the tasks execution (default OFF)" OFF)
if(VISUALIZE_EXECUTION)
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DUSE_VTK>)
  option(ENABLE_VELOCITY_RENDERING "Visualize the execution of velocity tasks (default OFF)" ON)
  option(ENABLE_STRESS_RENDERING "Visualize the execution of stress tasks (default OFF)" OFF)
  option(ENABLE_CLUSTER_RENDERING "Visualize the cluster activity during tasks execution (default OFF)" OFF)

  # Adjust compile options accordingly
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DOPENGL2_SUPPORTED>) ## Not working for Raspberry PI 3
  add_compile_options(-DVISUALIZE_EXECUTION)
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

option(DOWNLOAD_MISSING_DEPS "Enable automatic download of missing dependencies (default ON)" ON)
if (DOWNLOAD_MISSING_DEPS)
    # Download and install vcpkg (if not found) that will be used for installing dependencies
    set(VCPKG_ROOT ${CMAKE_SOURCE_DIR}/thirdparty/vcpkg/src)
    find_program(VCPKG_EXE vcpkg HINTS ${VCPKG_ROOT})
    if (NOT VCPKG_EXE)
        message(STATUS "vcpkg has not been found")
        message(STATUS "[START] vcpkg install")
        configure_file(${CMAKE_SOURCE_DIR}/cmake/resources/CMakeLists.txt.vcpkg ${VCPKG_ROOT}/../build/CMakeLists.txt)
        execute_process(COMMAND           ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
                        RESULT_VARIABLE   result
                        WORKING_DIRECTORY ${VCPKG_ROOT}/../build)
        if(result)
            message(FATAL_ERROR "CMake step for vcpkg failed: ${result}")
        endif()
        execute_process(COMMAND           ${CMAKE_COMMAND} --build .
                        RESULT_VARIABLE   result
                        WORKING_DIRECTORY ${VCPKG_ROOT}/../build)
        if(result)
            message(FATAL_ERROR "Build step for vcpkg failed: ${result}")
        endif()

        message(STATUS "[STOP] vcpkg install")
    else()
        get_filename_component(VCPKG_ROOT ${VCPKG_EXE} DIRECTORY)
        message(STATUS "Found vcpkg : ${VCPKG_EXE}")
    endif()
        
    set(CMAKE_TOOLCHAIN_FILE "${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake")
    message(STATUS "Using toolchain : ${CMAKE_TOOLCHAIN_FILE}")


    # set target platform
    if ("${CMAKE_HOST_SYSTEM_PROCESSOR}" MATCHES "(x86_64)|(AMD64)")
        set(HOST_ARCH "x64")
    elseif ("${CMAKE_HOST_SYSTEM_PROCESSOR}" MATCHES "(x86)|(i386)|(i686)")
        set(HOST_ARCH "x86")
    elseif ("${CMAKE_HOST_SYSTEM_PROCESSOR}" MATCHES "(aarch64)")
        set(HOST_ARCH "arm64")
    elseif ("${CMAKE_HOST_SYSTEM_PROCESSOR}" MATCHES "(arm)")
        set(HOST_ARCH "arm")
    else()
        message(FATAL_ERROR "Unsupported system : ${CMAKE_HOST_SYSTEM_PROCESSOR}")
    endif()

    string(TOLOWER ${CMAKE_SYSTEM_NAME} SYSTEM_NAME)
    set(VCPKG_TARGET_TRIPLET "${HOST_ARCH}-${SYSTEM_NAME}")

    # define helper macro for downloading packages through vpkg
    macro(_vcpkg_install package)
        message("[VCPKG] Installing package : ${package}:${TARGET_PLATFORM}")
        execute_process(COMMAND           vcpkg install ${package}:${TARGET_PLATFORM}
                        RESULT_VARIABLE   STATUS
                        WORKING_DIRECTORY ${VCPKG_ROOT})
        if (STATUS)
            message(FATAL_ERROR "Failed to install ${package}")
        else()
            message("[VCPKG] Successfully installed package : ${package}:${TARGET_PLATFORM}")
        endif()
    endmacro()
endif()

option(SEWAS_WITH_PARSEC "Use PaRSEC for task-scheduling (default ON)" ON)

if (SEWAS_WITH_PARSEC AND NOT SEWAS_DISTRIBUTED)
  message(FATAL_ERROR "The distributed version of SeWaS relies on PaRSEC.")
endif(SEWAS_WITH_PARSEC AND NOT SEWAS_DISTRIBUTED)

if (SEWAS_DISTRIBUTED)
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DSEWAS_DISTRIBUTED>)
endif(SEWAS_DISTRIBUTED)

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


# -------------------------------------------------------------------
# GTest
# -------------------------------------------------------------------
# Download and unpack googletest at configure time
if (BUILD_TESTING)
    configure_file(cmake/resources/CMakeLists.txt.gtest googletest-download/CMakeLists.txt)
    execute_process(COMMAND           ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
                    RESULT_VARIABLE   result
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
    if(result)
        message(FATAL_ERROR "CMake step for googletest failed: ${result}")
    endif()
    execute_process(COMMAND           ${CMAKE_COMMAND} --build .
                    RESULT_VARIABLE   result
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
    if(result)
        message(FATAL_ERROR "Build step for googletest failed: ${result}")
    endif()

    # Prevent overriding the parent project's compiler/linker
    # settings on Windows
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
    add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
    ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
    EXCLUDE_FROM_ALL)

    # The gtest/gtest_main targets carry header search path
    # dependencies automatically when using CMake 2.8.11 or
    # later. Otherwise we have to add them here ourselves.
    if (CMAKE_VERSION VERSION_LESS 2.8.11)
        include_directories("${gtest_SOURCE_DIR}/include")
    endif()
endif(BUILD_TESTING)

# -------------------------------------------------------------------
# spdlog
# -------------------------------------------------------------------
find_package(spdlog 1.3 CONFIG)
if (NOT spdlog_FOUND)
    # spdlog not found; trying to download it
    if (DOWNLOAD_MISSING_DEPS)
        _vcpkg_install(spdlog)
        find_package(spdlog 1.3 CONFIG REQUIRED)
    else()
        message(FATAL_ERROR "spdlog not found. Please set DOWNLOAD_MISSING_DEPS for enabling automatic download of missing dependencies")        
    endif()
endif()
set(SEWAS_LIBRARIES ${SEWAS_LIBRARIES} spdlog::spdlog)

# -------------------------------------------------------------------
# ADIOS2 is used for IO
# -------------------------------------------------------------------
if (ENABLE_IO)
  find_package(ADIOS2 REQUIRED)
  set(SEWAS_LIBRARIES ${SEWAS_LIBRARIES} ${ADIOS2_LIBRARIES})
endif(ENABLE_IO)

# -------------------------------------------------------------------
# MPI
# -------------------------------------------------------------------
if (SEWAS_DISTRIBUTED)
  # MPI
  find_package(MPI 2.0 REQUIRED)
  include_directories("${MPI_CXX_INCLUDE_PATH}")
  set(SEWAS_LIBRARIES ${SEWAS_LIBRARIES} ${MPI_CXX_LIBRARIES})
endif(SEWAS_DISTRIBUTED)

# -------------------------------------------------------------------
# PaRSEC
# -------------------------------------------------------------------
if(SEWAS_WITH_PARSEC)
  # PaRSEC
  find_package(PARSEC 2.0 REQUIRED)
  include_directories("${PARSEC_INCLUDE_DIRS}")
  link_directories("${PARSEC_LIBRARY_DIRS}")
  set(SEWAS_LIBRARIES ${SEWAS_LIBRARIES} ${PARSEC_LIBRARIES})
endif(SEWAS_WITH_PARSEC)

# -------------------------------------------------------------------
# Boost
# -------------------------------------------------------------------
find_package(Boost 1.62 COMPONENTS program_options)
if (NOT Boost_FOUND)
    # boost-program-options not found; trying to download it
    if (DOWNLOAD_MISSING_DEPS)
        _vcpkg_install(boost-program-options)
        _vcpkg_install(boost-property-tree)
        find_package(Boost 1.62 COMPONENTS program_options REQUIRED)
    else()
        message(FATAL_ERROR "Boost not found. Please set DOWNLOAD_MISSING_DEPS for enabling automatic download of missing dependencies")        
    endif()
endif()
include_directories(${Boost_INCLUDE_DIR})
set(SEWAS_LIBRARIES ${SEWAS_LIBRARIES} ${Boost_LIBRARIES})

# -------------------------------------------------------------------
# Eigen
# -------------------------------------------------------------------
find_package(Eigen3 3.3.3)
if (NOT Eigen3_FOUND)
    # eigen3 not found; trying to download it
    if (DOWNLOAD_MISSING_DEPS)
        _vcpkg_install(eigen3)
        find_package(Eigen3 3.3.3 REQUIRED)
    else()
        message(FATAL_ERROR "eigen3 not found. Please set DOWNLOAD_MISSING_DEPS for enabling automatic download of missing dependencies")        
    endif()
endif()
include_directories("${EIGEN3_INCLUDE_DIR}")

# -------------------------------------------------------------------
# Boost.SIMD
# -------------------------------------------------------------------
if (BOOST_SIMD_VECTORIZATION)
  find_package(Boost.SIMD 4.17 REQUIRED)
  include_directories("${Boost.SIMD_INCLUDE_DIRS}")
endif(BOOST_SIMD_VECTORIZATION)

# -------------------------------------------------------------------
# Intel TBB concurrent containers are used for collecting metrics
# -------------------------------------------------------------------
if (COLLECT_STATS)
    find_package(TBB 4.3)
    if (NOT TBB_FOUND)
        # TBB not found; trying to download it
        if (DOWNLOAD_MISSING_DEPS)
            _vcpkg_install(tbb)
            find_package(TBB 3.3.3 REQUIRED)
        else()
            message(FATAL_ERROR "TBB not found. Please set DOWNLOAD_MISSING_DEPS for enabling automatic download of missing dependencies")        
        endif()
    endif()
    include_directories("${TBB_INCLUDE_DIRS}")
    set(SEWAS_LIBRARIES ${SEWAS_LIBRARIES} ${TBB_LIBRARIES_RELEASE} ${TBB_LIBRARIES_DEBUG})
endif(COLLECT_STATS)

# -------------------------------------------------------------------
# matplotlibcpp
# -------------------------------------------------------------------
# The package matplotlib-cpp is used to plot data -> https://github.com/lava/matplotlib-cpp
# TODO Add an alternative for using xplot -> https://github.com/QuantStack/xplot
if (PLOT_DATA)
  include_directories("${CMAKE_SOURCE_DIR}/thirdparty/matplotlibcpp")
  find_package(PythonLibs 2.7)
  target_include_directories(sewas PRIVATE ${PYTHON_INCLUDE_DIRS})
  set(SEWAS_LIBRARIES ${SEWAS_LIBRARIES} ${PYTHON_LIBRARIES})
endif(PLOT_DATA)

# -------------------------------------------------------------------
# VTK is used to visualize the tasks execution
# -------------------------------------------------------------------
if (VISUALIZE_EXECUTION)
  find_package(VTK 8.1 REQUIRED)
  include_directories("${VTK_INCLUDE_DIRS}")
  set(SEWAS_LIBRARIES ${SEWAS_LIBRARIES} ${VTK_LIBRARIES})
endif(VISUALIZE_EXECUTION)

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
  src/MetricsManager.cxx
  src/ExecutionContext.cxx
  src/LogManager.cxx
  src/IOManager.cxx)

include_directories(include)

add_executable(sewas ${SOURCES} ${JDF2C_SOURCES})
set_property(TARGET sewas PROPERTY CXX_STANDARD 17)
target_link_libraries(sewas ${SEWAS_LIBRARIES})

# -------------------------------------------------------------------
# -------------------------------------------------------------------
if (BUILD_TESTING)
  # unit tests
  add_executable(FDOTest
    tests/FDOUnitTest.cxx
    src/CartesianMesh3D.cxx
    src/LogManager.cxx)
  add_executable(HaloExtractTest
    tests/HaloExtractUnitTest.cxx
    src/CartesianMesh3D.cxx
    src/Mesh3DPartitioning.cxx
    src/DataSet.cxx
    src/HaloManager.cxx
    src/LogManager.cxx)
  add_executable(HaloUpdateTest
    tests/HaloUpdateUnitTest.cxx
    src/CartesianMesh3D.cxx
    src/Mesh3DPartitioning.cxx
    src/DataSet.cxx
    src/HaloManager.cxx
    src/LogManager.cxx)
  add_executable(ComputeVelocityTest
    tests/ComputeVelocityUnitTest.cxx
    src/CartesianMesh3D.cxx
    src/Mesh3DPartitioning.cxx
    src/SEWASParameterManager.cxx
    src/SEWASPaRSEC.cxx
    src/ExecutionContext.cxx
    src/MetricsManager.cxx
    src/LinearSeismicWaveModel.cxx
    src/DataSet.cxx
    src/HaloManager.cxx
    src/LogManager.cxx
    src/IOManager.cxx
    ${JDF2C_SOURCES})

  set_property(TARGET FDOTest PROPERTY CXX_STANDARD 17)
  set_property(TARGET HaloExtractTest PROPERTY CXX_STANDARD 17)
  set_property(TARGET HaloUpdateTest PROPERTY CXX_STANDARD 17)
  set_property(TARGET ComputeVelocityTest PROPERTY CXX_STANDARD 17)

  target_link_libraries(FDOTest gtest_main ${SEWAS_LIBRARIES})
  target_link_libraries(HaloExtractTest gtest_main ${SEWAS_LIBRARIES})
  target_link_libraries(HaloUpdateTest gtest_main ${SEWAS_LIBRARIES})
  target_link_libraries(ComputeVelocityTest gtest_main ${SEWAS_LIBRARIES})
endif(BUILD_TESTING)
