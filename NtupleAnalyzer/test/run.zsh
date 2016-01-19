#!/bin/env zsh

cdir=$(pwd)/../
NtupleDir=$(pwd)/../../NtupleProducer/
export LD_LIBRARY_PATH=${cdir}:${NtupleDir}:${NtupleDir}/obj:$LD_LIBRARY_PATH

# tools: plot, tran

#./../Analyzer \
#--file input_WZJets_MC.txt \
#--tree Nt \
#--outfile ./output_WZJets_MC \
#--nmax -1 \
#--isdata 0 \
#--nowe 8235531 \
#--xsec 4.42965   \
#--lumi 1600 \

./../Analyzer \
--file input.txt \
--tree Nt \
--outfile ./output \
--nmax -1 \
--isdata 0 \
--nowe 1 \
--xsec 1 \
--lumi 1 \
