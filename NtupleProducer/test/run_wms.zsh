#!/bin/env zsh

anaName="ttH"
runName="selOpt"
voName="vo.sbg.in2p3.fr"

jname="${anaName}_${runName}"

mode="tzq"
nmax=-1

fpath="listsTZQ_temp/"

#if [[ $mode == "tth" ]]; then
#fpath="listsTTH/"
#elif [[ $mode == "tzq" ]]; then
#fpath="listsTZQ/"
#fi

fid="wmsid_${jname}.txt"
fre="wmsid_${jname}_RESUBMIT.txt"

doRe=0
if [[ -f ${fre} ]]; then
doRe=1
fi

rm -f ${fid}
rm -f logRF/*
wdir=SUBMITWMS
rm -rf ${wdir}
mkdir ${wdir}

#userareaSRM="srm://sbgse1.in2p3.fr:8446/dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/TOY_${jname}"
#userareaRF="/dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/TOY_${jname}"
#userareaRFLOG="/dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/LOG_${jname}"

#if [[ ${doRe} == 0 ]]; then
#/usr/bin/rfrm -rf ${userareaRF}
#/usr/bin/rfmkdir ${userareaRF}

#/usr/bin/rfrm -rf ${userareaRFLOG}
#/usr/bin/rfmkdir ${userareaRFLOG}
#fi

#userareaRF="${scrdir}/TOY_${jname}"
#userareaRFLOG="${scrdir}/LOG_${jname}"

#if [[ ${doRe} == 0 ]]; then
#rm -rf ${userareaRF}
#mkdir ${userareaRF}

#rm -rf ${userareaRFLOG}
#mkdir ${userareaRFLOG}
#fi

echo "prepare LIBs"
./prepareLib.zsh

glite-wms-job-delegate-proxy -e https://sbgwms2.in2p3.fr:7443/glite_wms_wmproxy_server -d kskovpen
#glite-wms-job-delegate-proxy -d kskovpen

stat="stat.txt"

#flist=$(ls ${fpath} | egrep 'NTuple_53X_DYJetsToLL_M-50_ID1')
flist=$(ls ${fpath})
farr=("${(@f)$(echo $flist)}")
fint=2
ffirst=1
echo $flist | while read line
do
sampleID=$(echo $line | sed 's%.txt%%g' | sed 's%_ID.*%%g')
sample=$(echo $line | sed 's%.txt%%g')
dataset=$(echo $sample | sed 's%_ID.*%%g')
fnext=$(echo ${farr[${fint}]} | sed 's%.txt%%g' | sed 's%_ID.*%%g')
statn=$(grep ${dataset} ${stat} | awk '{print $2}')

if [[ ${ffirst} == 1 ]]; then
mkdir ${wdir}/${sampleID}

#dire=$(rfdir ${userareaRF}/${sampleID} 2>&1)
direFAIL=$(echo $dire | grep "No such")
if [[ ${direFAIL} != "" && ${doRe} == 0 ]]; then
#/usr/bin/rfmkdir ${userareaRF}/${sampleID}
fi
jidx=0
ffirst=0
else
jidx=$[$jidx+1]
fi
fout=$(echo ${sampleID})
lout=$(echo ${sample})
cp j.jdl ${wdir}/${sampleID}/temp.jdl

cat ${wdir}/${sampleID}/temp.jdl | \
sed "s%line2%${line}%g" | \
sed "s%nmax%${nmax}%g" | \
sed "s%lout%${lout}%g" | \
sed "s%dataset%${dataset}%g" | \
sed "s%mode%${mode}%g" | \
sed "s%jname%${jname}%g" | \
sed "s%statn%${statn}%g" | \
sed "s%anaName%${anaName}%g" | \
sed "s%runName%${runName}%g" | \
sed "s%voName%${voName}%g" | \
sed "s%INPUTFILELIST%${fpath}${line}%g" \
> ${wdir}/${sampleID}/j_${jidx}.jdl

rm ${wdir}/${sampleID}/temp.jdl

if [[ ${fnext} != ${sampleID} ]]; then
rr=()
if [[ ${doRe} == 1 ]]; then
#rr=($(grep ${sampleID} ${fre} | awk '{print $2}'))
rr=($(grep ${sampleID} ${fre}))
fi
if [[ $#rr != 0 || ${doRe} == 0 ]]; then
echo "$sampleID"
#  mv ${wdir}/${sampleID} ${wdir}/${sampleID}TEMP
#  mkdir ${wdir}/${sampleID}
#  llist=($(ls ${wdir}/${sampleID}TEMP))
#    for si in $rr
#    do
#    rname="j_${si}.jdl"
#    sid=$[${si}+1]
#    rname=$llist[${sid}]
#    cp ${wdir}/${sampleID}TEMP/${rname} ${wdir}/${sampleID}/.
#    done
#  rm -rf ${wdir}/${sampleID}TEMP
#outp=$(glite-wms-job-submit -a --collection ${wdir}/${sampleID} 2>&1 | grep "https://sbgwms1.in2p3.fr:9000")
outp=$(glite-wms-job-submit -e https://sbgwms2.in2p3.fr:7443/glite_wms_wmproxy_server -a --collection ${wdir}/${sampleID} 2>&1 | grep "https://sbgwms2.in2p3.fr:9000")
#glite-wms-job-submit -a --collection ${wdir}/${sampleID}
if [[ $outp == "" ]]; then
echo "Input sandbox oversized"
exit 1
fi
echo "$sampleID   $outp" >> ${fid}
fi #rrdoRe
ffirst=1
fi

fint=$[$fint+1]

done
