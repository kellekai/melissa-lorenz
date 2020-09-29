#!/bin/bash
if [ -z "$MELISSA_DA_INSTALL" ]; then
  echo "please set MELISSA_DA_INSTALL first"
  exit 1
fi
make
cd melissa-pdaf
bash compile.sh
