cmake_minimum_required(VERSION 3.10.3)

project(parsec-download)

include(ExternalProject)

if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
    # Pinned to the last commit before PaRSEC's "Reshape of datacopies" rework
    # (8f813a179, first shipped in parsec-3.0.2012-rc1). That change and everything
    # built on top of it, all the way through parsec-4.0.2411, breaks cross-rank
    # halo exchange for SeWaS's task graph (MPI_ERR_TRUNCATE in the 2+ rank case;
    # confirmed by bisection, root cause not yet identified upstream).
    ExternalProject_Add(parsec
        GIT_REPOSITORY    https://github.com/icldisco/parsec
        GIT_TAG           4659c27061ffe2c23b034e09becf7e44e54ef6d8
        SOURCE_DIR        "${PARSEC_DOWNLOAD_ROOT}/src"
        BINARY_DIR        "${PARSEC_DOWNLOAD_ROOT}/build"
        PATCH_COMMAND     cd ${PARSEC_DOWNLOAD_ROOT}/src && git apply ${PARSEC_DOWNLOAD_ROOT}/patches/increase_max_inout_deps.patch && sed -i "s/parsec_profiling_stream_t \\*new )/parsec_profiling_stream_t *new_stream )/" parsec/profiling.h
        CONFIGURE_COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE} -DBUILD_TESTING=OFF -DPARSEC_DIST_WITH_MPI=${SEWAS_DISTRIBUTED} -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_C_FLAGS=-fcommon -DCMAKE_INSTALL_PREFIX=${THIRDPARTY_ROOT}/install ${PARSEC_DOWNLOAD_ROOT}/src
        BUILD_COMMAND     make -j
        INSTALL_COMMAND   make install
        TEST_COMMAND      ""
    )
else()
    message(FATAL_ERROR "PaRSEC based version of SeWaS is not yet available on this system: ${CMAKE_SYSTEM_NAME}")
endif()

