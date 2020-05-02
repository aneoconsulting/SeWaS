include(GoogleTest)

set_source_files_properties(${JDF2C_SOURCES}
                            PROPERTIES GENERATED TRUE)

### 1) Unit testing

add_executable(FDOTest
    FDOUnitTest.cxx
    ../src/LogManager.cxx)
add_executable(HaloExtractTest
    HaloExtractUnitTest.cxx
    ../src/CartesianMesh3D.cxx
    ../src/Mesh3DPartitioning.cxx
    ../src/SEWASParameterManager.cxx
    ../src/DataSet.cxx
    ../src/HaloManager.cxx
    ../src/LogManager.cxx)
add_executable(HaloUpdateTest
    HaloUpdateUnitTest.cxx
    ../src/CartesianMesh3D.cxx
    ../src/Mesh3DPartitioning.cxx
    ../src/SEWASParameterManager.cxx
    ../src/DataSet.cxx
    ../src/HaloManager.cxx
    ../src/LogManager.cxx)
add_executable(HaloExchangeTest
    HaloExchangeUnitTest.cxx
    ../src/CartesianMesh3D.cxx
    ../src/Mesh3DPartitioning.cxx
    ../src/SEWASParameterManager.cxx
    ../src/DataSet.cxx
    ../src/HaloManager.cxx
    ../src/LogManager.cxx)
add_executable(ComputeVelocityTest
    ComputeVelocityUnitTest.cxx
    ../src/CartesianMesh3D.cxx
    ../src/Mesh3DPartitioning.cxx
    ../src/SEWASParameterManager.cxx
    ../src/SEWASPaRSEC.cxx
    ../src/SEWASSequential.cxx
    ../src/ExecutionContext.cxx
    ../src/MetricsManager.cxx
    ../src/LinearSeismicWaveModel.cxx
    ../src/DataSet.cxx
    ../src/HaloManager.cxx
    ../src/LogManager.cxx
    ../src/IOManager.cxx
    ${JDF2C_SOURCES})

foreach(target FDOTest HaloExtractTest HaloUpdateTest HaloExchangeTest ComputeVelocityTest)
    set_property(TARGET ${target} PROPERTY CXX_STANDARD 17)
    target_link_libraries(${target} gtest_main ${SEWAS_LIBRARIES})

    # discovering tests
    #if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
    #    gtest_discover_tests(${target} EXTRA_ARGS "--gtest_output=xml:${CMAKE_CURRENT_BINARY_DIR}\\" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    #elseif ("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
    #    gtest_discover_tests(${target} EXTRA_ARGS "--gtest_output=xml:${CMAKE_CURRENT_BINARY_DIR}/" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    #else()
    #    message(WARNING "GTest output directory is not set for this system : ${CMAKE_SYSTEM_NAME}")
    #    gtest_discover_tests(${target} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    #endif()
    gtest_discover_tests(${target} EXTRA_ARGS "--gtest_output=xml")
endforeach()

# copy input data file used by some tests
configure_file(${CMAKE_SOURCE_DIR}/data/input/TestX.json ${CMAKE_CURRENT_BINARY_DIR}/TestX.json)
configure_file(${CMAKE_SOURCE_DIR}/data/input/TestC.json ${CMAKE_CURRENT_BINARY_DIR}/TestC.json)

### 2) Integration testing

add_test(NAME TestC
         COMMAND $<TARGET_FILE:sewas> --cx 51 --cy 37 --cz 23 --P 1 --Q 1 --R 1 --nthreads 2 --dfile=${CMAKE_SOURCE_DIR}/data/input/TestC.json)

find_package(Python3 COMPONENTS Interpreter)
if (Python3_FOUND)
    ## Check tiling (cx, cy, cz)
    add_test(NAME TestTiling
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/testers/test_tiling.py -v
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

    ## Check parallel partitionning (P, Q, R)
    if (SEWAS_DISTRIBUTED)
        add_test(NAME TestParallelPartitionning
                 COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/testers/test_parallel_partitionning.py -v
                 WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
        set_tests_properties(TestParallelPartitionning PROPERTIES ENVIRONMENT "MPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}")
        set_tests_properties(TestParallelPartitionning PROPERTIES ENVIRONMENT "MPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}")
    endif()

    if (ENABLE_IO)
        set_tests_properties(TestTiling PROPERTIES ENVIRONMENT "PYTHONPATH=${ADIOS2_DIR}/../../../bin")
    endif()
else()
    message(WARNING "Python3 has not been found. Integration testing will not run")
endif()
