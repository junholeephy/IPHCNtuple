#!/bin/bash

#eosdir="/store/user/chanon/TZQ/TestNtuplesV7_Syst"
dir=/tmp/chanon

for inputfile in `cat FileList_Syst.txt`
do

  echo ${dir}/${inputfile}
#  echo $line > tmp
#  proc=`awk '{print $1}' tmp`
#  inputfile=`awk '{print $2}' tmp`
#  echo $proc $inputfile

  root -l -q 'SplitSyst.C("'${dir}/${inputfile}'")' #| awk 'NR==3' > tmp
  #nEntries=`cat tmp`
  #echo nEntries=$nEntries

  #nEntriesTot=$(( nEntriesTot + nEntries ));
  #echo nEntriesTot=$nEntriesTot

done 

#cd /tmp/chanon/
#ls *.root > l

#for file in `cat l`
#do
#  cmsStage $file ${eosdir}/$file
#  echo $file
#done

#cd -

#echo nEntriesTot=$nEntriesTot
