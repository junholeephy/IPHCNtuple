#!/bin/sh

export X509_USER_PROXY=/home-pbs/nchanon/proxy/x509up_u8148

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /home-pbs/nchanon/CMSSW_7_4_12_patch4/src
export SCRAM_ARCH=slc6_amd64_gcc491
eval `scramv1 runtime -sh`
cd -

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${dout}/../:${dout}/../../NtupleProducer/:${dout}/../../../KinFit/


echo ${xsec}

line2=${line2}
fout=${fout}
isdata=${isdata}
nowe=${nowe}
xsec=${xsec}
dout=${dout}
dout_f=${dout_f}
sample=${sample}
lumi=${lumi}
dataset=${dataset}



echo "Executing .././NtupleAnalyzer --file ${line2} --outfile ${dout_f}${fout} --isdata ${isdata} --doSystCombine ${doSystCombine} --nowe ${nowe} --xsec ${xsec} --lumi ${lumi} --nmax ${nmax} --dataset ${dataset}"
${dout}/../Analyzer --file ${line2} --outfile ${dout_f}${fout} --isdata ${isdata} --doSystCombine ${doSystCombine} --nowe ${nowe} --xsec ${xsec} --lumi ${lumi} --nmax ${nmax} --tree Nt --dataset ${dataset}
