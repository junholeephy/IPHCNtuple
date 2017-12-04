#!/bin/sh

#export X509_USER_PROXY=/home-pbs/xcoubez/proxy/x509up_u7650
export X509_USER_PROXY=/home-pbs/ntonon/proxy/x509up_u8066

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /home-pbs/ntonon/tHq/CMSSW_8_0_20/src/
export SCRAM_ARCH=slc6_amd64_gcc530
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
