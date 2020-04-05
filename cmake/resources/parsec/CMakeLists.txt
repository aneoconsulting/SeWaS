cmake_minimum_required(VERSION 3.10.3)

project(parsec-download)

include(ExternalProject)

if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
    ExternalProject_Add(parsec
        GIT_REPOSITORY    https://bitbucket.org/icldistcomp/parsec.git
        GIT_TAG           8f72a0ecab1296fc3d90ba7c164695f4a3738950
        SOURCE_DIR        "${PARSEC_DOWNLOAD_ROOT}/src"
        BINARY_DIR        "${PARSEC_DOWNLOAD_ROOT}/build"
        PATCH_COMMAND     cd ${PARSEC_DOWNLOAD_ROOT}/src && git apply ${PARSEC_DOWNLOAD_ROOT}/patches/increase_max_inout_deps.patch
        CONFIGURE_COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE} -DBUILD_TESTING=OFF -DPARSEC_DIST_WITH_MPI=${SEWAS_DISTRIBUTED} -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_INSTALL_PREFIX=${THIRDPARTY_ROOT}/install ${PARSEC_DOWNLOAD_ROOT}/src
        BUILD_COMMAND     make -j
        INSTALL_COMMAND   make install
        TEST_COMMAND      ""
    )
else()
    message(FATAL_ERROR "PaRSEC based version of SeWaS is not yet available on this system: ${CMAKE_SYSTEM_NAME}")
endif()

