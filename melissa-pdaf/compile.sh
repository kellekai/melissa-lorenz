#!/bin/bash
#before using script export melissa install directory!
source $MELISSA_DA_INSTALL/bin/melissa-da_set_env.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MELISSA_DA_INSTALL/lib
export LIBRARY_PATH=$LIBRARY_PATH:$MELISSA_DA_INSTALL/lib
mkdir build; cd build
cmake \
  -DCMAKE_PREFIX_PATH=$MELISSA_DA_INSTALL \
  -DCMAKE_INSTALL_PREFIX=$(dirname `pwd`)/install \
  ..
make
make install
