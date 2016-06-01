#!/bin/bash

maindir=/opt/sbg/scratch1/cms/nchanon/ntuplesProd_763/toytoy_763_TTbar/

for dir in `cat ProcessList_TTbar.txt`
do

	ls ${maindir}/${dir} > ${dir}_tmp.txt
	for i in `cat ${dir}_tmp.txt`
	do
 		echo ${maindir}/${dir}/${i} 
	done > ${dir}.txt
	rm  ${dir}_tmp.txt

done

