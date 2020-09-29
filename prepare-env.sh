#!/bin/bash
if [ "${BASH_SOURCE[0]}" = "${0}" ]; then
  echo "this script needs to be sourced"
  exit 1
fi
if [ -z "$MELISSA_DA_INSTALL" ]; then
  echo "please set MELISSA_DA_INSTALL first"
  return
fi
source  $MELISSA_DA_INSTALL/bin/melissa-da_set_env.sh
export DATASET_PATH=`pwd`/data
