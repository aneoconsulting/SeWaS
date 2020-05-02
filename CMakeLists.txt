cmake_minimum_required(VERSION 3.10.3)

option(DOWNLOAD_MISSING_DEPS "Enable automatic download of missing dependencies (default ON)" ON)
if (DOWNLOAD_MISSING_DEPS)
    set(VCPKG_ROOT ${CMAKE_SOURCE_DIR}/thirdparty/vcpkg/src)

    find_program(VCPKG_EXE vcpkg HINTS ${VCPKG_ROOT})
    if (NOT VCPKG_EXE)
        message(FATAL_ERROR "DOWNLOAD_MISSING_DEPS is ON and we are unable to find vcpkg. Please download and install it or run bootstrap script before continuing.")
    else()
        get_filename_component(VCPKG_ROOT ${VCPKG_EXE} DIRECTORY)
        message(STATUS "Found vcpkg : ${VCPKG_EXE}")
    endif()
        
    # CMAKE_TOOLCHAIN_FILE must be set before the first call of project() directive
    # https://github.com/microsoft/vcpkg/blob/master/docs/users/integration.md#using-an-environment-variable-instead-of-a-command-line-option
    set(CMAKE_TOOLCHAIN_FILE ${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake CACHE STRING "Toolchain file used by vcpkg")
    message(STATUS "Using toolchain : ${CMAKE_TOOLCHAIN_FILE}")
endif()


project(SEWAS)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# version number
set(SEWAS_VERSION_MAJOR 0)
set(SEWAS_VERSION_MINOR 2)
set(SEWAS_VERSION_PATCH 0)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules" CACHE STRING "")


# ----------------------------------------------
# Build options
# ----------------------------------------------

option(SEWAS_DISTRIBUTED            "Build SeWaS for distributed memory platforms (default OFF)"                        OFF)
option(SEWAS_WITH_PARSEC            "Use PaRSEC for task-scheduling (default OFF)"                                      OFF)
option(VERBOSE                      "Enable collecting and displaying of traced events at runtime (default OFF)"        OFF)
option(BUILD_TESTING                "Build the testing tree (default OFF)"                                              OFF)
option(COLLECT_STATS                "Collect statistics during execution (default ON)"                                  ON)
option(SHOW_PER_THREAD_STATS        "Enable showing statistics at thread-level (default OFF)"                           OFF)
option(PLOT_DATA                    "Enable ploting the data of a tile object (default OFF)"                            OFF)
option(ENABLE_IO                    "Enable performing IO at runtime (default OFF)"                                     OFF)
option(VISUALIZE_EXECUTION          "Visualize the tasks execution (default OFF)"                                       OFF)
option(ENABLE_VELOCITY_RENDERING    "Visualize the execution of velocity tasks (default OFF)"                           ON)
option(ENABLE_STRESS_RENDERING      "Visualize the execution of stress tasks (default OFF)"                             OFF)
option(ENABLE_CLUSTER_RENDERING     "Visualize the cluster activity during tasks execution (default OFF)"               OFF)
option(USE_BLOCKWISE_FDO            "Use block operations for evaluating the finite difference operator (default ON)"   ON)
option(VECTORIZE_COMPUTATIONS       "Explicit vectorization of the computation (default ON)"                            ON)
option(EIGEN_VECTORIZATION          "Use Eigen as the vectorization engine (default ON)"                                ON)
option(BOOST_SIMD_VECTORIZATION     "Use Boost.SIMD as the vectorization engine (default OFF)"                          OFF)


# ----------------------------------------------
# Build configuration
# ----------------------------------------------

if (DOWNLOAD_MISSING_DEPS)
    set(THIRDPARTY_ROOT ${CMAKE_SOURCE_DIR}/thirdparty)

    # set vcpkg triplet for the target platform
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
    set(TARGET_PLATFORM "${HOST_ARCH}-${SYSTEM_NAME}")
    
    
    # define helper macro for downloading packages through vpckg
    macro(_vcpkg_install package)
        message("[VCPKG] Installing package : ${package}:${TARGET_PLATFORM}")
        execute_process(COMMAND           ${VCPKG_EXE} install ${package}:${TARGET_PLATFORM}
                        RESULT_VARIABLE   STATUS
                        WORKING_DIRECTORY ${VCPKG_ROOT})
        if (STATUS)
            message(FATAL_ERROR "Failed to install ${package} : ${STATUS}")
        else()
            message("[VCPKG] Successfully installed package : ${package}:${TARGET_PLATFORM}")
        endif()
    endmacro()
endif()

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

	# add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DEIGEN_ENABLE_AVX512>)
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND ${HOST_ARCH} MATCHES "(x64)")
    add_compile_options(-mcx16)
    add_compile_options(-m64)
endif()

if (SEWAS_DISTRIBUTED)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DSEWAS_DISTRIBUTED>)
endif(SEWAS_DISTRIBUTED)

if(SEWAS_WITH_PARSEC)
    add_compile_options(-DSEWAS_WITH_PARSEC)
endif()

if (VERBOSE)
    add_compile_options(-DVERBOSE)
endif(VERBOSE)

if (BUILD_TESTING)
    add_compile_options(-DBUILD_TESTING)
endif(BUILD_TESTING)

if (COLLECT_STATS)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DCOLLECT_STATS>)

    if(SHOW_PER_THREAD_STATS)
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DSHOW_PER_THREAD_STATS>)
    endif(SHOW_PER_THREAD_STATS)
endif(COLLECT_STATS)

if (PLOT_DATA)
    # Adjust compile options accordingly
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DUSE_MATPLOTLIB>)
endif(PLOT_DATA)

if (ENABLE_IO)
    # Adjust compile options accordingly
    add_compile_options(-DENABLE_IO)
endif(ENABLE_IO)

if(VISUALIZE_EXECUTION)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DUSE_VTK>)
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

if (USE_BLOCKWISE_FDO)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DBLOCKWISE_FDO>)
endif(USE_BLOCKWISE_FDO)

if(VECTORIZE_COMPUTATIONS)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DVECTORIZE_COMPUTATIONS>)

    if (EIGEN_VECTORIZATION)
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DEIGEN_VECTORIZATION>)
    endif(EIGEN_VECTORIZATION)

    if (BOOST_SIMD_VECTORIZATION)
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-DBOOST_SIMD_VECTORIZATION>)
    endif(BOOST_SIMD_VECTORIZATION)
endif(VECTORIZE_COMPUTATIONS)


# ----------------------------------------------
# Build configuration completed
# ----------------------------------------------

# -------------------------------------------------------------------
# GTest
# -------------------------------------------------------------------
# Download and unpack googletest at configure time
if (BUILD_TESTING)
    configure_file(cmake/resources/gtest/CMakeLists.txt googletest-download/CMakeLists.txt)
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
set(SEWAS_LIBRARIES ${SEWAS_LIBRARIES} spdlog::spdlog spdlog::spdlog_header_only)

# -------------------------------------------------------------------
# MPI
# -------------------------------------------------------------------
if (SEWAS_DISTRIBUTED)
    # MPI
    find_package(MPI 3.0)
    if (NOT MPI_FOUND)
        # openmpi not found; trying to download it
        if (DOWNLOAD_MISSING_DEPS)
            _vcpkg_install(openmpi)
            find_package(MPI 3.0 REQUIRED)
        else()
            message(FATAL_ERROR "MPI not found. Please set DOWNLOAD_MISSING_DEPS for enabling automatic download of missing dependencies")        
        endif()
    endif()
    include_directories("${MPI_CXX_INCLUDE_PATH}")
    set(SEWAS_LIBRARIES ${SEWAS_LIBRARIES} ${MPI_CXX_LIBRARIES})
endif(SEWAS_DISTRIBUTED)

# -------------------------------------------------------------------
# ADIOS2 is used for IO
# -------------------------------------------------------------------
if (ENABLE_IO)
    find_package(ADIOS2 PATHS ${THIRDPARTY_ROOT}/install/lib/cmake/adios2)
    if (NOT ADIOS2_FOUND)
        message(WARNING "ADIOS2 not found")
        if (DOWNLOAD_MISSING_DEPS)
            # ADIOS2 not found; trying to download it
            set(ADIOS_DOWNLOAD_ROOT ${THIRDPARTY_ROOT}/adios)
            configure_file(${CMAKE_SOURCE_DIR}/cmake/resources/adios/CMakeLists.txt ${ADIOS_DOWNLOAD_ROOT}/build/adios-download/CMakeLists.txt)
            execute_process(COMMAND           ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
                            RESULT_VARIABLE   result
                            WORKING_DIRECTORY ${ADIOS_DOWNLOAD_ROOT}/build/adios-download)
            if(result)
                message(FATAL_ERROR "CMake step for ADIOS2 failed: ${result}")
            endif()
            execute_process(COMMAND           ${CMAKE_COMMAND} --build .
                            RESULT_VARIABLE   result
                            WORKING_DIRECTORY ${ADIOS_DOWNLOAD_ROOT}/build/adios-download)
            if(result)
                message(FATAL_ERROR "Build step for ADIOS2 failed: ${result}")
            endif()
                        
            find_package(ADIOS2 REQUIRED PATHS ${THIRDPARTY_ROOT}/install/lib/cmake/adios2)
        else(DOWNLOAD_MISSING_DEPS)
            message(FATAL_ERROR "ADIOS2 not found. Please set DOWNLOAD_MISSING_DEPS for enabling automatic download of missing dependencies")
        endif(DOWNLOAD_MISSING_DEPS)
    else(NOT ADIOS2_FOUND)
        message(STATUS "ADIOS2 found")
        message(STATUS "  ADIOS2_DIR = ${ADIOS2_DIR}")
    endif(NOT ADIOS2_FOUND)
    set(SEWAS_LIBRARIES ${SEWAS_LIBRARIES} ${ADIOS2_LIBRARIES})
    
    if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
        add_compile_definitions(NOMINMAX)
    endif()
endif(ENABLE_IO)

# -------------------------------------------------------------------
# PaRSEC
# -------------------------------------------------------------------
if(SEWAS_WITH_PARSEC)
    # PaRSEC
    if (NOT SEWAS_DISTRIBUTED)
        message(FATAL_ERROR "The distributed version of SeWaS relies on PaRSEC.")
    endif(NOT SEWAS_DISTRIBUTED)

    find_package(PaRSEC CONFIG PATHS ${THIRDPARTY_ROOT}/install/share/cmake/parsec)
    if (NOT PARSEC_DIR)
        message(WARNING "PaRSEC not found")
        if (DOWNLOAD_MISSING_DEPS)
            # PaRSEC not found; trying to download it
            set(PARSEC_DOWNLOAD_ROOT ${THIRDPARTY_ROOT}/parsec)
            configure_file(${CMAKE_SOURCE_DIR}/cmake/resources/parsec/CMakeLists.txt ${PARSEC_DOWNLOAD_ROOT}/build/parsec-download/CMakeLists.txt)
            execute_process(COMMAND           ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
                            RESULT_VARIABLE   result
                            WORKING_DIRECTORY ${PARSEC_DOWNLOAD_ROOT}/build/parsec-download)
            if(result)
                message(FATAL_ERROR "CMake step for parsec failed: ${result}")
            endif()
            execute_process(COMMAND           ${CMAKE_COMMAND} --build .
                            RESULT_VARIABLE   result
                            WORKING_DIRECTORY ${PARSEC_DOWNLOAD_ROOT}/build/parsec-download)
            if(result)
                message(FATAL_ERROR "Build step for parsec failed: ${result}")
            endif()
                        
            find_package(PaRSEC CONFIG REQUIRED PATHS ${THIRDPARTY_ROOT}/install/share/cmake/parsec)
        else(DOWNLOAD_MISSING_DEPS)
            message(FATAL_ERROR "PaRSEC not found. Please set DOWNLOAD_MISSING_DEPS for enabling automatic download of missing dependencies")
        endif(DOWNLOAD_MISSING_DEPS)
    endif(NOT PARSEC_DIR)


    # Compilation rules of JDF files
    include(RulesJDF)
    find_program(parsec-ptgpp_EXE parsec-ptgpp PATHS ${PARSEC_BINARY_DIRS})
    if (NOT parsec-ptgpp_EXE)
        message(FATAL_ERROR "parsec-ptgpp not found")
    else()
        set(PARSEC_DEP_MGT "HASH_TABLE" CACHE STRING "PaRSEC dependency management (INDEX_ARRAY|HASH_TABLE)")
        mark_as_advanced(FORCE PARSEC_DEP_MGT)
        if ("${PARSEC_DEP_MGT}" STREQUAL "INDEX_ARRAY")
            set(PARSEC_PTGPP_CFLAGS "-M;index-array") # Workaround suggested by George to avoid sporadic non-termination
        endif()

        jdf_rules(JDF2C_SOURCES ${CMAKE_SOURCE_DIR}/src/runtimes/parsec/sewas.jdf)
        include_directories(${CMAKE_CURRENT_BINARY_DIR})
    endif()

    
    message(STATUS "PaRSEC found")
    message(STATUS "  PARSEC_DIR = ${PARSEC_DIR}")
    message(STATUS "  PARSEC_INCLUDE_DIRS = ${PARSEC_INCLUDE_DIRS}")
    message(STATUS "  PARSEC_LIBRARIES = ${PARSEC_LIBRARIES}")
    message(STATUS "  parsec-ptgpp_EXE = ${parsec-ptgpp_EXE}")

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
    message(DEPRECATION "Boot.SIMD backend is not fully tested. Its support will be discarded in a future release")
    find_package(Boost.SIMD 4.17 REQUIRED)
    include_directories("${Boost.SIMD_INCLUDE_DIRS}")
endif(BOOST_SIMD_VECTORIZATION)

# -------------------------------------------------------------------
# matplotlibcpp
# -------------------------------------------------------------------
# The package matplotlib-cpp is used to plot data -> https://github.com/lava/matplotlib-cpp
# TODO Add an alternative for using xplot -> https://github.com/QuantStack/xplot
if (PLOT_DATA)
    find_package(Python2 COMPONENTS Development)
    if (NOT Python2_FOUND)
        # PythonLibs not found; trying to download it
        if (DOWNLOAD_MISSING_DEPS)
            _vcpkg_install(python2)
            _vcpkg_install(matplotlib-cpp)
            find_package(Python2 COMPONENTS Development REQUIRED)
        else()
            message(FATAL_ERROR "Python2 not found. Please set DOWNLOAD_MISSING_DEPS for enabling automatic download of missing dependencies")        
        endif()
    endif()
    set(SEWAS_INCLUDE_DIRS ${SEWAS_INCLUDE_DIRS} ${Python2_INCLUDE_DIRS})
    set(SEWAS_LIBRARIES ${SEWAS_LIBRARIES} Python2::Python)
    add_compile_definitions(WITHOUT_NUMPY)
endif(PLOT_DATA)

# -------------------------------------------------------------------
# VTK is used to visualize the tasks execution
# -------------------------------------------------------------------
if (VISUALIZE_EXECUTION)
    find_package(VTK 8.1)
    if (NOT VTK_FOUND)
        # VTK not found; trying to download it
        if (DOWNLOAD_MISSING_DEPS)
            _vcpkg_install(vtk)
            find_package(VTK 8.1 CONFIG REQUIRED)
        else()
            message(FATAL_ERROR "VTK not found. Please set DOWNLOAD_MISSING_DEPS for enabling automatic download of missing dependencies")        
        endif()
    endif()
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
target_include_directories(sewas PRIVATE ${SEWAS_INCLUDE_DIRS})
target_link_libraries(sewas ${SEWAS_LIBRARIES})

# -------------------------------------------------------------------
# -------------------------------------------------------------------
if (BUILD_TESTING)
    include(CTest)

    # unit tests
    add_subdirectory(tests)
endif(BUILD_TESTING)
