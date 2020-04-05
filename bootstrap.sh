#!/bin/sh -e

SEWAS_ROOT=$(X= cd -- "$(dirname -- "$0")" && pwd -P)
. $SEWAS_ROOT/scripts/bootstrap-sewas.sh $SEWAS_ROOT
