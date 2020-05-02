cmake_minimum_required(VERSION 3.10.3)

project(adios-download)

include(ExternalProject)

ExternalProject_Add(adios
    GIT_REPOSITORY    https://github.com/ornladios/ADIOS2.git
    GIT_TAG           v2.5.0
    SOURCE_DIR        "${ADIOS_DOWNLOAD_ROOT}/src"
    BINARY_DIR        "${ADIOS_DOWNLOAD_ROOT}/build"
    PATCH_COMMAND     ""
    CONFIGURE_COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE} -DADIOS2_BUILD_EXAMPLES=OFF -DADIOS2_BUILD_TESTING=OFF -DADIOS2_USE_MPI=${SEWAS_DISTRIBUTED} -DADIOS2_USE_Python=ON -DCMAKE_INSTALL_PREFIX=${THIRDPARTY_ROOT}/install ${ADIOS_DOWNLOAD_ROOT}/src
    BUILD_COMMAND     make -j
    INSTALL_COMMAND   make install
    TEST_COMMAND      ""
)

