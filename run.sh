#!/bin/bash

export OMP_NUM_THREADS=8

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
${DIR}/src/reference
