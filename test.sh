#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd ${DIR}

# Compile testdriver
cd testdriver
make
cd ..

# Compile program
./compile.sh

# Run test driver with basic test case
testdriver/testdriver ./run.sh basic.test
