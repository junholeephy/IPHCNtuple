#!/bin/sh

export X509_USER_PROXY=/home-pbs/lebihan/someone/proxy/x509up_u6155

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /home-pbs/lebihan/someone/CMSSW_7_4_12_patch4/src
export SCRAM_ARCH=slc6_amd64_gcc491
eval `scramv1 runtime -sh`
cd -

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${dout}/../

line2=${line2}
fout=${fout}
isdata=${isdata}
noe=${noe}
xsec=${xsec}
dout=${dout}
dout_f=${dout_f}
sample=${sample}


echo "Executing .././NtupleProducer --file ${line2} --outfile ${dout_f}${fout} --isdata ${isdata} --noe ${noe} --xsec ${xsec} --nmax ${nmax}"
${dout}/./NtupleProducer --file ${line2} --outfile ${dout_f}${fout} --isdata ${isdata} --noe ${noe} --xsec ${xsec} --nmax ${nmax} --tree FlatTree/tree
