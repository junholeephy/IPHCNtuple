#!/bin/env zsh

# don't forget /opt/sbg/scratch1/cms


echo "Don't forget to update the lumi and the maximum number of events to run on in this script if needed !"

isdata=0
lumi=2260

cp /tmp/x509up_u8148 /home-pbs/nchanon/proxy/.

jName=${1}
if [[ ${jName} == "" ]]; then
  echo "Please specify the run name"
  exit 1
fi

que="sbg_local"

export HOME=$(pwd)

dout="/home-pbs/nchanon/CMSSW_7_4_12_patch4/src/IPHCNtuple/NtupleAnalyzer/test/"
dout_f="/home-pbs/nchanon/CMSSW_7_4_12_patch4/src/IPHCNtuple/NtupleAnalyzer/test/"

runName="toy_${jName}"
logName="log_${jName}"

rm -rf ${runName}
mkdir ${runName}
rm -rf ${logName}
mkdir ${logName}

nmax=-1

fxsec="table_MC_Akoula-patch2_02022016.txt"

fdir=$(ls -d lists*)

echo $fdir

echo $fdir | while read line
do
fpath="${HOME}/${line}/"
flist=$(ls ${fpath})
dir=${line}

echo $flist | while read line
do
  jidx=0
  sample=$(echo $line | sed 's%.txt%%g')
  dataset=$(echo $sample | sed 's%_ID..*%%g')
  if [[ ! -d ${runName}/${dataset} ]]; then
    mkdir ${runName}/${dataset}
  fi
  linexsec=$(grep $dataset $fxsec)
  nowe=$(echo $linexsec | awk '{print $3}')
  xsec=$(echo $linexsec | awk '{printf $2}')
  if [[ $nowe == "" ]]; then
    nowe=1
  fi
  if [[ $xsec == "" ]]; then
    xsec=1
  fi
  fl=$(echo $sample | cut -c1-1)
  fl2=$(echo $sample | cut -c1-11)
  datamc=""
  if [[ $fl == "J" ]]; then
    isdata=1
    datamc="DATA"
    nmax=${nmax}
  else
    isdata=0
    datamc="MC"
    nmax=${nmax}
  fi
  
  
  fout=`echo ${runName}/${dataset}/${line}_${jidx} | sed 's%.txt%%g'`
  lout=$(echo ${line}_${jidx} | sed 's%.txt%%g')

  echo "${dataset}: $nowe $xsec $lumi"
  #echo "${fpath}${line}"
 
  qsub -N ${dir} -q ${que} -o ${logName}/${sample}.log -j oe single_batch_job.sh \
-v dout=${dout},line2=${fpath}${line},dout_f=$dout_f$,fout=${fout},nowe=${nowe},xsec=${xsec},lumi=${lumi},isdata=${isdata},dataset=${dataset},nmax=${nmax}
done

echo "going to sleep 2700 s (45 mn)"
sleep 2700

done 
