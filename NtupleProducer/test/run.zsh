#!/bin/env zsh

cdir=$(pwd)/../
export LD_LIBRARY_PATH=${cdir}:${cdir}/obj:$LD_LIBRARY_PATH

infl="input.txt"

./NtupleProducer     \
--file ${infl}       \
--tree FlatTree/tree \
--nmax 1000
