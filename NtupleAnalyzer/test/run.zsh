#!/bin/env zsh

cdir=$(pwd)/../
NtupleDir=$(pwd)/../../NtupleProducer/
KinFitDir=$(pwd)/../../../KinFit/
export LD_LIBRARY_PATH=${cdir}:${NtupleDir}:${KinFitDir}:${NtupleDir}/obj:$LD_LIBRARY_PATH

# tools: plot, tran

./../Analyzer \
--file test.txt \
--tree Nt \
--outfile ./output_ttH_MC \
--nmax 10000 \
--isdata 0 \
--doSystCombine 1 \
--nowe 1 \
--xsec 1   \
--lumi 1 \

#./../Analyzer \
#--file backup_lists_akoula_patch5_20161003_MC/ttHToNonbb_M125_13TeV_powheg_pythia8.txt \
#--tree Nt \
#--outfile ./output_ttH_MC \
#--nmax 10000 \
#--isdata 0 \
#--nowe 1 \
#--xsec 1   \
#--lumi 1 \

#./../Analyzer \
#--file input.txt \
#--tree Nt \
#--outfile ./output \
#--nmax -1 \
#--isdata 0 \
#--nowe 1 \
#--xsec 1 \
#--lumi 1 \
