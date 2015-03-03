#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd ${DIR}/src

#chmod +x ./lib/tbb/bin/tbbvars.sh
#./lib/tbb/bin/tbbvars.sh intel64

make clean
make
