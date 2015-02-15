#!/bin/bash
# Make script executable from any location
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd ${DIR}

# Compile logic
tar -czvf submission.tar.gz compile.sh run.sh src