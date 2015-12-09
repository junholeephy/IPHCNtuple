#!/bin/bash

opt=x3_LargeRange_FullStat

for i in `seq 0 179`
do
  if [ ! -e Jobs_${opt}/output_TTW_x3_LargeRange_FullStat_${i}.root ] ; then echo bsub -q 2nd RunBatchMEM_TTW_${i} ; fi
done

for i in `seq 0 165`
do
  if [ ! -e Jobs_${opt}/output_TTZ_x3_LargeRange_FullStat_${i}.root ] ; then echo bsub -q 2nd RunBatchMEM_TTZandGstar_${i} ; fi
done

for i in `seq 0 221`
do
  if [ ! -e Jobs_${opt}/output_TTG_x3_LargeRange_FullStat_${i}.root ] ; then echo bsub -q 2nd RunBatchMEM_TTG_${i} ; fi
done

for i in `seq 0 166`
do
  if [ ! -e Jobs_${opt}/output_TT_x3_LargeRange_FullStat_${i}.root ] ; then echo bsub -q 2nd RunBatchMEM_TT_${i} ; fi
done

for i in `seq 0 75`
do
  if [ ! -e Jobs_${opt}/output_WZ_x3_LargeRange_FullStat_${i}.root ] ; then echo bsub -q 2nd RunBatchMEM_WZ_${i} ; fi
done

for i in `seq 0 1880`
do
  if [ ! -e Jobs_${opt}/output_TTH_x3_LargeRange_FullStat_${i}.root ] ; then echo bsub -q 2nd RunBatchMEM_TTHflsl_${i} ; fi
done
