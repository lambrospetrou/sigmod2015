#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

export GOPATH=${DIR}

cd ${DIR}/src/validator
go clean
go build -v