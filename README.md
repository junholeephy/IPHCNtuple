ttH analysis framework
============

based on IPHCFlatTree

Install
-------

```
# ROOT is required

git clone https://github.com/kskovpen/ttH.git

cd ttH/

NtupleProducer/ - analysis package to produce user class-based ntuples with all object
parameters needed to perform ttH selection

cd NtupleProducer/
make
cd test/
./run.zsh

NtupleAnalyzer/ - analysis package to apply ttH selection based on
NtupleProducer ntuples

cd NtupleAnalyzer/
make
cd test/
./run.zsh

```
