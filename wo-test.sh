#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd ${DIR}

# Compile testdriver
#echo -e "Compiling testdriver"
#cd testdriver
#make
#cd ..

# Compile program
#echo -e "Compiling program"
#./compile.sh

# Run test driver with basic test case
#testdriver/testdriver ./run.sh small.test
echo -e "Running program"
testdriver/testdriver ./run.sh $1
