#!/bin/bash

# HTSlib first
cd htslib
if [ ! -f libhts.a ]; then
  make
fi
cd $DEPS
