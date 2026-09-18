cmake_minimum_required(VERSION 3.10.3)

project(starpu-download)

include(ExternalProject)

if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
    # SEWAS_STARPU_ENABLE_FXT/FXT_ROOT are substituted to literal values (e.g.
    # "if(ON)", a literal prefix path) by the configure_file() call in the
    # top-level CMakeLists.txt that generates this file from the template.
    set(STARPU_EXTRA_CONFIGURE_ARGS "")
    if(${SEWAS_STARPU_ENABLE_FXT})
        list(APPEND STARPU_EXTRA_CONFIGURE_ARGS "--with-fxt=${FXT_ROOT}")
    endif()

    # Release tarball (ships a pre-generated `configure`, so no autoreconf/automake/
    # libtool is required at build time). Version/hash pinned to StarPU 1.4.7.
    #
    # --enable-maxcpus=128: StarPU's default (`--enable-maxcpus=auto`) bakes in
    # the CPU count of the machine running THIS configure step, rounded up to
    # the next power of 2 -- not the machine(s) StarPU will actually run on.
    # Confirmed on the graviton_tests cluster: StarPU built on a 4-vCPU
    # instance capped STARPU_MAXCPUS at 4, so a StarPU built there silently
    # ran only 4 CPU workers even on an 8-vCPU node ("7 CPU cores devices
    # requested. Only 4 enabled." at startup). Pin it explicitly so the
    # binary isn't limited by whatever machine happened to build it.
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
                          --enable-maxcpus=128
                          ${STARPU_EXTRA_CONFIGURE_ARGS}
        BUILD_COMMAND     make -j
        INSTALL_COMMAND   make install
        TEST_COMMAND      ""
    )
else()
    message(FATAL_ERROR "StarPU based version of SeWaS is not yet available on this system: ${CMAKE_SYSTEM_NAME}")
endif()
