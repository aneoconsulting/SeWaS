cmake_minimum_required(VERSION 3.10.3)

project(vcpkg-download)

include(ExternalProject)

if ("${CMAKE_HOST_SYSTEM_PROCESSOR}" STREQUAL "aarch64")
    set(_VCPKG_BOOTSTRAP_OPTIONS "-useSystemBinaries")
endif()

if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
    ExternalProject_Add(vcpkg
        GIT_REPOSITORY    https://github.com/microsoft/vcpkg.git
        GIT_TAG           40f3c709db80acf154ac4b17a1f83c564ebd022e
        SOURCE_DIR        "../src"
        BUILD_IN_SOURCE   true
        CONFIGURE_COMMAND ""
        BUILD_COMMAND     .\\bootstrap-vcpkg.bat
        INSTALL_COMMAND   ""
        TEST_COMMAND      ""
    )
elseif ("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
    # arm64-linux support is native to vcpkg since well before this pin; the old
    # add_aarch64_support.patch/community triplet copy is both unneeded and no
    # longer applies (scripts/buildsystems/vcpkg.cmake has since been rewritten
    # upstream).
    ExternalProject_Add(vcpkg
        GIT_REPOSITORY    https://github.com/microsoft/vcpkg.git
        GIT_TAG           40f3c709db80acf154ac4b17a1f83c564ebd022e
        SOURCE_DIR        "../src"
        BUILD_IN_SOURCE   true
        CONFIGURE_COMMAND ""
        BUILD_COMMAND     ./bootstrap-vcpkg.sh ${_VCPKG_BOOTSTRAP_OPTIONS}
        INSTALL_COMMAND   ""
        TEST_COMMAND      ""
    )
else()
    message(FATAL_ERROR "Unsupported system: ${CMAKE_SYSTEM_NAME}")
endif()

