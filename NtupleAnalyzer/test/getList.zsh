#!/bin/env zsh

dirName=${1}

#if [[ ${dirName} == "" ]]; then
#  echo "Specify toy directory"
#  exit 1
#fi

#fpath="/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_12/src/tHFCNC/NtupleProducer/test/toy${dirName}/"
fpath="/opt/sbg/scratch1/cms/xcoubez/ntuples_prod_medusa_patch4/toyMC_${dirName}/"

liDATA=($(ls ${fpath} | grep "Run2016"))
liMC=($(ls ${fpath}))

nFilesDATA=2 # per job
nFilesMC=1
outDir="lists"

rm -rf ${outDir}
mkdir ${outDir}

rm -f /tmp/tempDATA.txt
for line in $liDATA
do
d1=$(echo $line)
liDATA2=$(ls ${fpath}${d1})
echo $liDATA2 | while read line2
do
f1=$(echo $line2)
file=$(echo ${fpath}${d1}/${f1})
echo "${file}" >> /tmp/tempDATA.txt
done
split -a 5 -l ${nFilesDATA} -d /tmp/tempDATA.txt /tmp/${d1}_
lsfi=$(ls /tmp/${d1}_*)
jid=0
echo $lsfi | while read fil
do
mv ${fil} ${outDir}/${d1}_ID${jid}.txt
jid=$[$jid+1]
done
rm -f /tmp/tempDATA.txt
done

rm -f /tmp/tempMC.txt
for line in $liMC
do
d1=$(echo $line)
liMC2=$(ls ${fpath}${d1})
echo $liMC2 | while read line2
do
f1=$(echo $line2)
file=$(echo ${fpath}${d1}/${f1})
echo "${file}" >> /tmp/tempMC.txt
done
nFiles=${nFilesMC}
split -a 5 -l ${nFiles} -d /tmp/tempMC.txt /tmp/${d1}_
lsfi=$(ls /tmp/${d1}_*)
jid=0
echo $lsfi | while read fil
do
mv ${fil} ${outDir}/${d1}_ID${jid}.txt
jid=$[$jid+1]
done
rm -f /tmp/tempMC.txt
done

