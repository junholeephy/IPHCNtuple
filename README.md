# IPHCNtuple
Ntuple production and analysis code

Install
-------

```
# The code is based on ROOT6, setup it manually or via CMSSW:

mkdir MyAnalysis/; cd MyAnalysis/
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_8_0_20
cd CMSSW_8_0_20/src/
cmsenv

# get the code from GIT
git clone -b Moriond2017  https://github.com/IPHC/ttH


cd IPHCNtuple/

# NtupleProducer: produce Ntuples from FlatTrees

cd NtupleProducer/
make
cd ../../

# NtupleAnalyzer: create histograms, TTrees, ASCII, etc output from Ntuples

cd ttH/NtupleAnalyzer/
make


# to commit & push 
git push origin Moriond2017
