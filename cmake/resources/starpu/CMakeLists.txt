cmake_minimum_required(VERSION 3.10.3)

project(starpu-download)

include(ExternalProject)

if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
    # Release tarball (ships a pre-generated `configure`, so no autoreconf/automake/
    # libtool is required at build time). Version/hash pinned to StarPU 1.4.7.
    ExternalProject_Add(starpu
        URL               https://files.inria.fr/starpu/starpu-1.4.7/starpu-1.4.7.tar.gz
        URL_HASH          SHA256=1eb3df551089153fe6e0b172ad95118434b4a81ea9eaa5a2c38725d0db53b13e
        SOURCE_DIR        "${STARPU_DOWNLOAD_ROOT}/src"
        BUILD_IN_SOURCE   1
        CONFIGURE_COMMAND <SOURCE_DIR>/configure
                          --prefix=${THIRDPARTY_ROOT}/install
                          --with-mpicc=${MPI_C_COMPILER}
                          --with-mpicxx=${MPI_CXX_COMPILER}
                          --disable-build-tests
                          --disable-build-examples
                          --disable-build-doc
                          --disable-fortran
                          --disable-cuda
                          --disable-opencl
        BUILD_COMMAND     make -j
        INSTALL_COMMAND   make install
        TEST_COMMAND      ""
    )
else()
    message(FATAL_ERROR "StarPU based version of SeWaS is not yet available on this system: ${CMAKE_SYSTEM_NAME}")
endif()
