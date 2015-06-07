#!/bin/env zsh

#export LD_PRELOAD=/usr/lib64/libglobus_gssapi_gsi.so.4

cdir=$(pwd)/../
export LD_LIBRARY_PATH=${cdir}:$LD_LIBRARY_PATH

#infl="lists/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_ID0.txt"
infl="input.txt"

./NtupleProducer \
--file ${infl} \
--tree FlatTree/tree
