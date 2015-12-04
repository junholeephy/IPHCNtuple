#!/bin/env zsh

cdir=$(pwd)/../
NtupleDir=$(pwd)/../../NtupleProducer/
export LD_LIBRARY_PATH=${cdir}:${NtupleDir}:${NtupleDir}/obj:$LD_LIBRARY_PATH

# tools: plot, tran

./../Analyzer \
--file input.txt \
--tree Nt \
--outfile ./ \
--nmax -1 \
--isdata 0 \
--nowe -1 \
--xsec -1. \
--lumi 2110 \
