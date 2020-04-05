cmake_minimum_required(VERSION 3.10.3)

project(vcpkg-download)

include(ExternalProject)

if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
    ExternalProject_Add(vcpkg
        GIT_REPOSITORY    https://github.com/microsoft/vcpkg.git
        GIT_TAG           master
        SOURCE_DIR        "../src"
        BUILD_IN_SOURCE   true
        CONFIGURE_COMMAND ""
        BUILD_COMMAND     .\\bootstrap-vcpkg.bat
        INSTALL_COMMAND   ""
        TEST_COMMAND      ""
    )
elseif ("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
    ExternalProject_Add(vcpkg
        GIT_REPOSITORY    https://github.com/microsoft/vcpkg.git
        GIT_TAG           master
        SOURCE_DIR        "../src"
        BUILD_IN_SOURCE   true
        CONFIGURE_COMMAND ""
        BUILD_COMMAND     ./bootstrap-vcpkg.sh
        INSTALL_COMMAND   ""
        TEST_COMMAND      ""
    )
else()
    message(FATAL_ERROR "Unsupported system: ${CMAKE_SYSTEM_NAME}")
endif()

