# NtProd - NtAnalyzer Workflow

README for the tHq-2016 branch, describing the basic steps to run the production chain from FlatTrees to ntuples ready for MVA analysis.

*Do not forget to source :*
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /cvmfs/cms.cern.ch/crab3/crab.sh
```

## FlatTreeProducer

(( Follow instructions from [IPHCFlatTree's README](https://github.com/IPHC/IPHCFlatTree/tree/tHq) ))


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
cd /home-pbs/username/MyAnalysis/CMSSW_8_0_20/src/ttH/NtupleProducer/src
```

* **NtupleProducer.cxx** - modify path to JEC files for data and MC, e.g. : 

```
jesTotal = new JetCorrectionUncertainty(*(new JetCorrectorParameters("/home-pbs/ntonon/tHq/CMSSW_8_0_25/src/IPHCFlatTree/FlatTreeProducer/data/jecFiles/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_UncertaintySources_AK4PFchs.txt", "Total")));
```


```
cd /home-pbs/username/MyAnalysis/CMSSW_8_0_20/src/ttH/NtupleProducer/test
```

* **input.txt** - include the the Flat Tree(s) path(s), e.g. : 

```
root://sbgse1.in2p3.fr//dpm/in2p3.fr/home/cms/phedex/store/user/XXX/output_*.root
```


* **input.txt** - include the the Flat Tree(s) path(s), e.g. : 

```
root://sbgse1.in2p3.fr//dpm/in2p3.fr/home/cms/phedex/store/user/XXX/output_*.root
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

//Or directly (e.g.) : 
./NtupleProducer --file input.txt  --outfile output --tree FlatTree/tree --nmax -1 --isdata 0
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

