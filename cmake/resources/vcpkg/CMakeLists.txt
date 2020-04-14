cmake_minimum_required(VERSION 3.10.3)

project(vcpkg-download)

include(ExternalProject)

if ("${CMAKE_HOST_SYSTEM_PROCESSOR}" STREQUAL "aarch64")
    set(_VCPKG_BOOTSTRAP_OPTIONS "-useSystemBinaries")
endif()

if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
    ExternalProject_Add(vcpkg
        GIT_REPOSITORY    https://github.com/microsoft/vcpkg.git
        GIT_TAG           941d5464544eb5812b2cbd6e687c3ebc78ed2624
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
        GIT_TAG           941d5464544eb5812b2cbd6e687c3ebc78ed2624
        SOURCE_DIR        "../src"
        BUILD_IN_SOURCE   true
        CONFIGURE_COMMAND ""
        PATCH_COMMAND     cp ../triplets/arm64-linux.cmake ../src/triplets/community && cd ../src && git apply ../patches/add_aarch64_support.patch
        BUILD_COMMAND     ./bootstrap-vcpkg.sh ${_VCPKG_BOOTSTRAP_OPTIONS}
        INSTALL_COMMAND   ""
        TEST_COMMAND      ""
    )
else()
    message(FATAL_ERROR "Unsupported system: ${CMAKE_SYSTEM_NAME}")
endif()

