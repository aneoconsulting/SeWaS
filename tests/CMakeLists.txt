include(GoogleTest)

add_executable(FDOTest
    FDOUnitTest.cxx
    ../src/CartesianMesh3D.cxx
    ../src/LogManager.cxx)
add_executable(HaloExtractTest
    HaloExtractUnitTest.cxx
    ../src/CartesianMesh3D.cxx
    ../src/Mesh3DPartitioning.cxx
    ../src/DataSet.cxx
    ../src/HaloManager.cxx
    ../src/LogManager.cxx)
add_executable(HaloUpdateTest
    HaloUpdateUnitTest.cxx
    ../src/CartesianMesh3D.cxx
    ../src/Mesh3DPartitioning.cxx
    ../src/DataSet.cxx
    ../src/HaloManager.cxx
    ../src/LogManager.cxx)
add_executable(HaloExchangeTest
    HaloExchangeUnitTest.cxx
    ../src/CartesianMesh3D.cxx
    ../src/Mesh3DPartitioning.cxx
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
    if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
        gtest_discover_tests(${target} EXTRA_ARGS "--gtest_output=xml:${CMAKE_CURRENT_BINARY_DIR}\\" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    elseif ("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
        gtest_discover_tests(${target} EXTRA_ARGS "--gtest_output=xml:${CMAKE_CURRENT_BINARY_DIR}/" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    else()
        message(WARNING "GTest output directory is not set for this system : ${CMAKE_SYSTEM_NAME}")
        gtest_discover_tests(${target} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    endif()
endforeach()

# add hard link to input data file used by some tests
# file(CREATE_LINK ${CMAKE_SOURCE_DIR}/data/input/TestX.json ${CMAKE_CURRENT_BINARY_DIR}/TestX.json COPY_ON_ERROR)
