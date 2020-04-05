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
endforeach()
