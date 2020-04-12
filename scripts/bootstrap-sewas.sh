#!/bin/sh -e

SEWAS_ROOT=$1

VCPKG_ROOT=$SEWAS_ROOT/thirdparty/vcpkg/src
VCPKG_BUILD=$SEWAS_ROOT/thirdparty/vcpkg/build

echo "[START] Bootstraping"

_vcpkg=`command -v $VCPKG_ROOT/vcpkg`
if [[ -x "$VCPKG_ROOT/vcpkg" ]]; then
	echo "vcpkg has not been found"
	
	echo "[START] install vcpkg"

	rm -rf $VCPKG_BUILD && mkdir -p $VCPKG_BUILD

	cp $SEWAS_ROOT/cmake/resources/vcpkg/CMakeLists.txt $VCPKG_BUILD

	cd $VCPKG_BUILD && pwd && cmake . && cmake --build . && cd -

	echo "[STOP] install vcpkg"
else
	echo "vcpkg found : " $_vcpkg
	$_vcpkg version
	$_vcpkg update
fi

echo "[STOP] Bootstraping"
