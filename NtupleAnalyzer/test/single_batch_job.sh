#!/bin/sh

export X509_USER_PROXY=/home-pbs/lebihan/someone/proxy/x509up_u6155

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /home-pbs/lebihan/someone/CMSSW_8_0_9/src
export SCRAM_ARCH=slc6_amd64_gcc530
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
