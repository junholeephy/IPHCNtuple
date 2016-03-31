#!/bin/bash

while read line
do

  echo $line > tmp
  proc=`awk '{print $1}' tmp`
  inputfile=`awk '{print $2}' tmp`
  echo $proc $inputfile

  root -l -q 'ReadEntries.C("'${inputfile}'")' | awk 'NR==3' > tmp
  nEntries=`cat tmp`
echo nEntries=$nEntries

done < FileList.txt
