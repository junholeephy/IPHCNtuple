# IPHCNtuple
Ntuple production and analysis code

Install
-------

```
# The code is based on ROOT6, setup it manually or via CMSSW:

mkdir MyAnalysis/; cd MyAnalysis/
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_8_0_9
cd CMSSW_8_0_9/src/
cmsenv

# get the code from GIT
git clone -b 80X-branch  https://github.com/IPHC/ttH


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

cd ttH/NtupleAnalyzer/
make


# to commit & push 
git push origin 80X-branch
