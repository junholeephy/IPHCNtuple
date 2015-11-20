# IPHCNtuple
Ntuple production and analysis code

Install
-------

```
# get the code from GIT
git clone https://github.com/IPHC/IPHCNtuple

cd IPHCNtuple/

# The code is based on ROOT6, setup it manually or via CMSSW:

mkdir MyAnalysis/; cd MyAnalysis/
export SCRAM_ARCH=slc6_amd64_gcc491
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_7_4_12_patch4
cd CMSSW_7_4_12_patch4/src/
cmsenv

# NtupleProducer: produce Ntuples from FlatTrees

cd NtupleProducer/
make
cd test/
./run.zsh

# NtupleAnalyzer: create histograms, TTrees, ASCII, etc output from Ntuples

cd NtupleAnalyzer/
make
cd test/
./run.zsh
