# IPHCNtuple
Ntuple production and analysis code

Install
-------

```
# The code is based on ROOT6, setup it manually or via CMSSW:

mkdir MyAnalysis/; cd MyAnalysis/
export SCRAM_ARCH=slc6_amd64_gcc493
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_7_6_3
cd CMSSW_7_6_3/src/
cmsenv

# get the code from GIT
git clone https://github.com/IPHC/IPHCNtuple

cd IPHCNtuple/

# NtupleProducer: produce Ntuples from FlatTrees

cd NtupleProducer/
make
cd ../../

# NtupleAnalyzer: create histograms, TTrees, ASCII, etc output from Ntuples

git clone https://github.com/kskovpen/KinFit
cd KinFit/
make
cd ../

cd IPHCNtuple/NtupleAnalyzer/
make
