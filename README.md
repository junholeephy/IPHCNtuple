# FlatTreeProducer - NtProd - NtAnalyzer Workflow

README for the tHq-2016 branch, describing the basic steps to run the production chain from datasets to ntuples ready for MVA analysis.

*Do not forget to source :*
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /cvmfs/cms.cern.ch/crab3/crab.sh
```

## FlatTreeProducer

Installing and running the IPHCFlatTree code to produce Flat Trees. Using tag "Walrus-patch-2".

### Installation

```
cd /home-pbs/username/
mkdir MyAnalysis
cd MyAnalysis

# CMSSW Release
RELEASE=8_0_25

# Setup release
cmsrel CMSSW_$RELEASE
cd CMSSW_X_Y_Z/src
cmsenv
git cms-init

# Clone this repo
git clone https://github.com/IPHC/IPHCFlatTree.git
cd IPHCFlatTree
git tag -l //List the available tags
git checkout tags/Walrus-patch2 //Use Walrus-patch2 tag instead of master branch
cd ..

# Egamma
git cms-merge-topic shervin86:Moriond2017_JEC_energyScales
cd EgammaAnalysis/ElectronTools/data; git clone https://github.com/ECALELFS/ScalesSmearings; cd -
git cms-merge-topic ikrav:egm_id_80X_v2

# Add MET filters
git cms-merge-topic -u cms-met:fromCMSSW_8_0_20_postICHEPfilter

# Tools needed for AK10 jet collection
git clone https://github.com/cms-jet/JetToolbox JMEAnalysis/JetToolbox 

# Add DeepCSV tagger
git cms-merge-topic -u mverzett:DeepFlavour-from-CMSSW_8_0_21
mkdir RecoBTag/DeepFlavour/data/; cd RecoBTag/DeepFlavour/data/; wget http://home.fnal.gov/~verzetti//DeepFlavour/training/DeepFlavourNoSL.json; cd -

# Compile the monster
scram b -j5
```

(( Instructions taken from [IPHCFlatTree's README](https://github.com/IPHC/IPHCFlatTree/tree/Walrus-patch2) ))


### Set-up


```
cd IPHCFlatTree/FlatTreeProducer/test/PROD
```
* **list.txt** - create it and add all the dataset names, e.g. : 

```
...
/THQ_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
...
```

* **submit.zsh** - modify the following :

```
...
slist="list.txt"
ver="name_of_your_ouput_dir"
...
```

* **crabConfigTemplate.py** - modify the following :

```
...
isData=0 #For MC
...
config.Data.splitting = 'FileBased' #For MC
#config.Data.splitting = 'LumiBased' #For data
...
```


### Launch the jobs

```
./submit.zsh
```


## IPHCNtuple

Ntuple Production & Analysis codes

### Installation

```
cd /home-pbs/username
mkdir MyAnalysis/; cd MyAnalysis/

#Could install it in 8_0_25, as IPHCFlatTree ?
cmsrel CMSSW_8_0_20
cd CMSSW_8_0_20/src/
cmsenv

# get the code from GIT
git clone -b tHq2016  https://github.com/IPHC/ttH


cd IPHCNtuple/

# NtupleProducer: produce Ntuples from FlatTrees

cd NtupleProducer/
make
cd ../../

# NtupleAnalyzer: create histograms, TTrees, ASCII, etc output from Ntuples

cd ttH/NtupleAnalyzer/
make


# to commit & push 
git push origin tHq2016
```

### Set-up & Run NtupleProducer


```
cd /home-pbs/username/MyAnalysis/CMSSW_8_0_20/src/ttH/NtupleProducer/test
```

* **input.txt** - include the the Flat Tree(s) path(s), e.g. : 

```
root://sbgse1.in2p3.fr//dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/FlatTree/Walrus-patch2-v20170615/THQ_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM/170615_174912/0000/output_1.root
```

* **split_into_lists.zsh** - modify the following lines with your own proxy, username, directory :

```
...
export x509_USER_PROXY=/home-pbs/ntonon/proxy/x509up_u8066
...
fpath="/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/FlatTree/output_dir/"
...
```


* **single_batch_job.zsh** - idem : 

```
...
export X509_USER_PROXY=/home-pbs/ntonon/proxy/x509up_u8066
...
cd /home-pbs/ntonon/tHq/CMSSW_8_0_20/src/
...
```

* **run_batch.zsh** - idem : 

```
...
cp /tmp/x509up_u8066 /home-pbs/ntonon/proxy
...
dout="/home-pbs/ntonon/tHq/CMSSW_8_0_20/src/ttH/NtupleProducer/test"
dout_f="/opt/sbg/scratch1/cms/ntonon/ntuples_prod_walrus_patch2/"
...
```

### Run

*Interactively*

```
./run.zsh
```

*Launch jobs*

```
./run_batch.zsh
```


### Set-up & Run NtupleAnalyzer


```
cd /home-pbs/username/MyAnalysis/CMSSW_8_0_20/src/ttH/NtupleAnalyzer/test
```
* **table.txt** - add list of samples, with cross section (pb-1) and number of weights (nowe), e.g. : 

```
...
THQ_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM_0000       0.7927       3495652
...
```


* **split_into_lists.zsh** - add path of directory containing input files, from NtupleProducer : 

```
...
fpath="/opt/sbg/scratch1/cms/ntonon/ntuples_prod_walrus_patch2/toyallStat/"
...
```

* **run_batch.zsh** - update lumi, proxy, username, directories : 

```
...
lumi=35900
...
cp /tmp/x509up_u8066 /home-pbs/ntonon/proxy
...
dout="/home-pbs/ntonon/tHq/CMSSW_8_0_20/src/ttH/NtupleAnalyzer/test"
dout_f="/opt/sbg/scratch1/cms/ntonon/Analyzer_ntuples_prod_walrus_patch2/"
...
fdir=$(ls -d lists_NameOfYourList)
...
```

-- For interactive running : 

* **run.zsh** - update lumi, xsec, nof weights, ... 

* **input.txt** - add path of input files, from NtupleProducer (for interactive running) : 


### Run

*Interactively*

```
./run.zsh
```

*Launch jobs*

```
./run_batch.zsh
```

