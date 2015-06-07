#!/bin/sh

export LFC_HOST=sbglfc1.in2p3.fr

lfile="tthrun.log"

echo "=====  Begin  =====" | tee ${lfile}
date | tee -a ${lfile}
echo "The program is running on $HOSTNAME" | tee -a ${lfile}
date | tee -a ${lfile}
echo "=====  End  =====" | tee -a ${lfile}

lsb_release -a | tee -a ${lfile}

WDIR=$(pwd)

line2=${1}
nmax=${2}
dataset=${3}
lout=${4}
mode=${5}
jname=${6}
statn=${7}
anaName=${8}
runName=${9}
voName=${10}

odir="/opt/sbg/scratch1/${voName}/kskovpen/${anaName}/${runName}"

export ROOTSYS=/cvmfs/cms.cern.ch/slc6_amd64_gcc472/lcg/root/5.32.00
ls $ROOTSYS/bin/thisroot.sh
source $ROOTSYS/bin/thisroot.sh
rootV=$(root-config --version)
echo "ROOT v${rootV} has been set up"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${WDIR}

#export LD_LIBRARY_PATH=\
#/cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/gcc/4.7.2/lib64:\
#/usr/lib64:\
#/cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/gcc/4.7.2/lib:\
#/cvmfs/cms.cern.ch/slc6_amd64_gcc472/cms/cmssw-patch/CMSSW_5_3_16_patch1/external/slc6_amd64_gcc472/lib/:\
#/cvmfs/cms.cern.ch/slc6_amd64_gcc472/cms/cmssw-patch/CMSSW_5_3_16_patch1/lib/slc6_amd64_gcc472/:\
#$LD_LIBRARY_PATH

export LD_LIBRARY_PATH=\
/cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/gcc/4.7.2/lib:\
/cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/gcc/4.7.2/lib64:\
/cvmfs/cms.cern.ch/slc6_amd64_gcc472/cms/cmssw-patch/CMSSW_5_3_16_patch1/external/slc6_amd64_gcc472/lib/:\
/cvmfs/cms.cern.ch/slc6_amd64_gcc472/cms/cmssw-patch/CMSSW_5_3_16_patch1/lib/slc6_amd64_gcc472/:\
/cvmfs/cms.cern.ch/slc6_amd64_gcc472/cms/cmssw/CMSSW_5_3_16/lib/slc6_amd64_gcc472/:\
$LD_LIBRARY_PATH

export PATH=/cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/gcc/4.7.2/bin:$PATH

# this allows rfio access on sl6
export LD_PRELOAD=/usr/lib64/libglobus_gssapi_gsi.so.4

chmod a+x tthrun

#tar -zxvf PtRelFall12.tar.gz
tar -zxvf LIB.tar.gz
mv TEMPLIB/* .
rm -rf TEMPLIB
tar -zxvf INPUT.tar.gz
mv TEMPINPUT/* .
rm -rf TEMPINPUT

mkdir xml
mv TTbarHAnalysis.xml xml/
mv TZqAnalysis.xml xml/

ls

##ldd tthrun | tee -a ${lfile}
##ldd libNTuple.so | tee -a ${lfile}

#/usr/bin/./make 2>&1 | tee -a ${lfile}
#chmod a+x tthrun

echo "Executing ../tthrun ${line2} ${nmax} ${dataset} ${mode} ${statn}" | tee -a ${lfile}
./tthrun ${line2} ${nmax} ${dataset} ${mode} ${statn} 2>&1 | tee -a ${lfile}
#/usr/bin/./make
#chmod a+x readFile
#fname=$(cat $line2) 
#ldd readFile | tee -a ${lfile} 
#./readFile "${fname}" 2>&1 | tee -a ${lfile}

#srmcp -retry_num 4 -retry_timeout 30000 file:///${WDIR}/output.root srm://sbgse1.in2p3.fr:8446/dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/TOY_${jname}/${dataset}/$lout.root
cp ${WDIR}/output.root ${odir}DATA/${dataset}/$lout.root

#if [ $? -ne 0 ]; then
#echo "retrying srmcp pour root"
#srmrm srm://sbgse1.in2p3.fr:8446/dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/TOY_${jname}/${dataset}/$lout.root
#srmcp -retry_num 4 -retry_timeout 30000 file:///${WDIR}/output.root srm://sbgse1.in2p3.fr:8446/dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/TOY_${jname}/${dataset}/$lout.root
#fi

#srmcp -retry_num 4 -retry_timeout 30000 file:///${WDIR}/${lfile} srm://sbgse1.in2p3.fr:8446/dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/LOG_${jname}/${lout}.log
cp ${WDIR}/${lfile} ${odir}LOG/${lout}.log

#if [ $? -ne 0 ]; then
#echo "retrying srmcp pour log"
#srmrm srm://sbgse1.in2p3.fr:8446/dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/LOG_${jname}/${lout}.log
#srmcp -retry_num 4 -retry_timeout 30000 file:///${WDIR}/${lfile} srm://sbgse1.in2p3.fr:8446/dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/LOG_${jname}/${lout}.log
#fi

rm -f ${lfile}
rm -f output.root
rm -rf electronSF.root
rm -rf muonSF.root
rm pileup_hadd_official_upto173692.root
rm tthrun
rm $line2
rm LIB.tar.gz
rm INPUT.tar.gz
rm *.so

