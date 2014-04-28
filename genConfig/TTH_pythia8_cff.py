import FWCore.ParameterSet.Config as cms

source = cms.Source("EmptySource")

generator = cms.EDFilter("Pythia8GeneratorFilter",
pythiaPylistVerbosity = cms.untracked.int32(0),
filterEfficiency = cms.untracked.double(1.0),
pythiaHepMCVerbosity = cms.untracked.bool(False),
comEnergy = cms.double(8000.0),
maxEventsToPrint = cms.untracked.int32(0),
crossSection  = cms.untracked.double(2.072e-03),
PythiaParameters = cms.PSet(
processParameters = cms.vstring(
'HiggsSM:gg2Httbar = on',
'HiggsSM:qqbar2Httbar = on',
'25:m0 = 125.0', # MadGraph
#'25:mWidth = 0.006382339', # MadGraph
#'25:doForceWidth = true',
'25:onMode = off',
'25:onIfAny = 24 -24'
#        '24:onMode = off',
#        '24:onIfAny = -11 -13 -15 11 13 15',
#'6:m0 = 173.0', # MadGraph
#'6:mWidth = 1.491500', # MadGraph
#'6:doForceWidth = true',
#'6:onMode = off',                                                                                                  
#'6:onIfAny = 24 -24'
),
parameterSets = cms.vstring(
'processParameters')
)
)

configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    name = cms.untracked.string('$Source: ,v $'),
    annotation = cms.untracked.string('PYTHIA8 TTH at sqrt(s) = 8TeV')
    )

ProductionFilterSequence = cms.Sequence(generator)


