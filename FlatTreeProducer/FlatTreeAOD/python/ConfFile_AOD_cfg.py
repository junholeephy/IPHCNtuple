import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
import os

options = VarParsing('analysis')
options.register('isData',False,VarParsing.multiplicity.singleton,VarParsing.varType.int,'Run on real data')
process = cms.Process("FlatTree")
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'root://sbgse1.in2p3.fr//dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/testFiles/MiniAOD/ttH_ev.root'
#        'root://sbgse1.in2p3.fr//dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/testFiles/TZqSynchMoreStat/001C0B68-536A-E311-B25F-002590D0B066.root'
    )
)

process.FlatTree = cms.EDAnalyzer('FlatTreeAOD'
)


process.TFileService = cms.Service("TFileService",
fileName = cms.string("output.root")
#fileName = cms.string(options.outputFile)
)

process.p = cms.Path(
                     process.FlatTree)
