#!/bin/bash

for opt in x3_LargeRange_FullStat
do

  for i in `seq 0 179` # 147 885
  do
    i1=$((i*6))
    i2=$((i*6+6))
    cp RunBatchMEM_template Jobs_${opt}/RunBatchMEM_TTW_$i 
    sed -i s/NMIN/$i1/g Jobs_${opt}/RunBatchMEM_TTW_$i
    sed -i s/NMAX/$i2/g Jobs_${opt}/RunBatchMEM_TTW_$i
    sed -i s/NUM/$i/g Jobs_${opt}/RunBatchMEM_TTW_$i
    sed -i s/OPTION/${opt}/g Jobs_${opt}/RunBatchMEM_TTW_$i
    sed -i s/PROC/TTW/g Jobs_${opt}/RunBatchMEM_TTW_$i
    echo bsub -q 2nd RunBatchMEM_TTW_$i
  done

  for i in `seq 0 165` #142 856
  do
    i1=$((i*6))
    i2=$((i*6+6))
    cp RunBatchMEM_template Jobs_${opt}/RunBatchMEM_TTZandGstar_$i        
    sed -i s/NMIN/$i1/g Jobs_${opt}/RunBatchMEM_TTZandGstar_$i
    sed -i s/NMAX/$i2/g Jobs_${opt}/RunBatchMEM_TTZandGstar_$i
    sed -i s/NUM/$i/g Jobs_${opt}/RunBatchMEM_TTZandGstar_$i
    sed -i s/OPTION/${opt}/g Jobs_${opt}/RunBatchMEM_TTZandGstar_$i
    sed -i s/PROC/TTZ/g Jobs_${opt}/RunBatchMEM_TTZandGstar_$i
    echo bsub -q 2nd RunBatchMEM_TTZandGstar_$i
  done

  for i in `seq 0 221` 
  do
    i1=$((i*6))
    i2=$((i*6+6))
    cp RunBatchMEM_template Jobs_${opt}/RunBatchMEM_TTG_$i
    sed -i s/NMIN/$i1/g Jobs_${opt}/RunBatchMEM_TTG_$i
    sed -i s/NMAX/$i2/g Jobs_${opt}/RunBatchMEM_TTG_$i
    sed -i s/NUM/$i/g Jobs_${opt}/RunBatchMEM_TTG_$i
    sed -i s/OPTION/${opt}/g Jobs_${opt}/RunBatchMEM_TTG_$i
    sed -i s/PROC/TTG/g Jobs_${opt}/RunBatchMEM_TTG_$i
    echo bsub -q 2nd RunBatchMEM_TTG_$i
  done

  for i in `seq 0 166` 
  do
    i1=$((i*6))
    i2=$((i*6+6))
    cp RunBatchMEM_template Jobs_${opt}/RunBatchMEM_TT_$i
    sed -i s/NMIN/$i1/g Jobs_${opt}/RunBatchMEM_TT_$i
    sed -i s/NMAX/$i2/g Jobs_${opt}/RunBatchMEM_TT_$i
    sed -i s/NUM/$i/g Jobs_${opt}/RunBatchMEM_TT_$i
    sed -i s/OPTION/${opt}/g Jobs_${opt}/RunBatchMEM_TT_$i
    sed -i s/PROC/TT/g Jobs_${opt}/RunBatchMEM_TT_$i
    echo bsub -q 2nd RunBatchMEM_TT_$i
  done

  for i in `seq 0 75`
  do
    i1=$((i*6))
    i2=$((i*6+6))
    cp RunBatchMEM_template Jobs_${opt}/RunBatchMEM_WZ_$i
    sed -i s/NMIN/$i1/g Jobs_${opt}/RunBatchMEM_WZ_$i
    sed -i s/NMAX/$i2/g Jobs_${opt}/RunBatchMEM_WZ_$i
    sed -i s/NUM/$i/g Jobs_${opt}/RunBatchMEM_WZ_$i
    sed -i s/OPTION/${opt}/g Jobs_${opt}/RunBatchMEM_WZ_$i
    sed -i s/PROC/WZ/g Jobs_${opt}/RunBatchMEM_WZ_$i
    echo bsub -q 2nd RunBatchMEM_TT_$i
  done

  for i in `seq 0 1880` #`seq 0 333`
  do
    i1=$((i*6)) #20, 6
    i2=$((i*6+6))
    cp RunBatchMEM_template Jobs_${opt}/RunBatchMEM_TTHflsl_${i}
    sed -i s/NMIN/$i1/g Jobs_${opt}/RunBatchMEM_TTHflsl_$i
    sed -i s/NMAX/$i2/g Jobs_${opt}/RunBatchMEM_TTHflsl_$i
    sed -i s/NUM/$i/g Jobs_${opt}/RunBatchMEM_TTHflsl_$i
    sed -i s/OPTION/${opt}/g Jobs_${opt}/RunBatchMEM_TTHflsl_$i
    sed -i s/PROC/TTH/g Jobs_${opt}/RunBatchMEM_TTHflsl_$i
    echo bsub -q 2nd RunBatchMEM_TTHflsl_$i
  done

done 
