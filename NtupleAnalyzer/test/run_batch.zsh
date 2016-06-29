#!/bin/env zsh

# don't forget /opt/sbg/scratch1/cms


echo "Don't forget to update the lumi and the maximum number of events to run on in this script if needed !"

isdata=0
doSystCombine=0
lumi=2070

cp /tmp/x509up_u6155 /home-pbs/lebihan/someone/proxy/.

jName=${1}
if [[ ${jName} == "" ]]; then
  echo "Please specify the run name"
  exit 1
fi

que="cms"

export HOME=$(pwd)

dout="/home-pbs/lebihan/someone/medusa_patch1_prod/ttH/NtupleAnalyzer/test/"
dout_f="/opt/sbg/scratch1/cms/lebihan/trees_analyzer_prod_medusa_patch1_v3/"

runName="toy_${jName}"
logName="log_${jName}"

rm -rf ${dout_f}/${runName}
mkdir ${dout_f}/${runName}
rm -rf ${logName}
mkdir ${logName}

nmax=-1

fxsec="table_MC_Medusa-patch1_20160624.txt"

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
    mkdir ${dout_f}/${runName}/${dataset}
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
    doSystCombine=0
    datamc="DATA"
    nmax=${nmax}
  else
    isdata=0
    doSystCombine=0
    datamc="MC"
    nmax=${nmax}
  fi
  
  
  fout=`echo ${runName}/${dataset}/${line}_${jidx} | sed 's%.txt%%g'`
  lout=`echo ${line}_${jidx} | sed 's%.txt%%g'`
  #fout=$(echo ${runName}/${dataset}/${line}_${jidx} | sed 's%.txt%%g')
  #lout=$(echo ${line}_${jidx} | sed 's%.txt%%g')

  echo "${dataset}: $nowe $xsec $lumi"
  #echo "${fpath}${line}"
 
  qsub -N ${dir} -q ${que} -o ${logName}/${sample}.log -j oe single_batch_job.sh \
-v dout=${dout},line2=${fpath}${line},dout_f=${dout_f},fout=${fout},nowe=${nowe},xsec=${xsec},lumi=${lumi},isdata=${isdata},doSystCombine=${doSystCombine},dataset=${dataset},nmax=${nmax}
done

echo "going to sleep 2700 s (45 mn)"
sleep 2700

done 
