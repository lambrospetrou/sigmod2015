#!/bin/bash

#export OMP_NUM_THREADS=`grep siblings < /proc/cpuinfo | head -1 | sed 's/.*\([0-9]\+\)/\1/'`

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
${DIR}/src/reference
