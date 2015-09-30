import FWCore.ParameterSet.Config as cms


process = cms.Process("ICALIB")
process.source = cms.Source("EmptyIOVSource",
    firstValue = cms.uint64(245807),#243664),#(243257), #(239640),
    lastValue = cms.uint64(245807),#243664), #(243257), #(239640),
    timetype = cms.string('runnumber'),
    interval = cms.uint64(1)
)



process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.threshold = cms.untracked.string('DEBUG')
process.MessageLogger.cout.threshold = cms.untracked.string('DEBUG')
process.MessageLogger.debugModules = cms.untracked.vstring("*")
process.MessageLogger.destinations = cms.untracked.vstring('cout')
process.MessageLogger.cout = cms.untracked.PSet( threshold = cms.untracked.string('DEBUG'))

#to configure conditions via the Global Tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'GR_E_V43', '') #''GR_P_V50', '')

#process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')

#Setup the SiSTripFedCabling and the SiStripDetCabling
#process.load("CondCore.DBCommon.CondDBCommon_cfi")
#process.CondDBCommon.connect='frontier://FrontierProd/CMS_COND_31X_STRIP'

#process.poolDBESSource = cms.ESSource( 'PoolDBESSource',
#                                       process.CondDBCommon,
#                                       BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
#                                       toGet = cms.VPSet( 
#                                                          cms.PSet( record = cms.string('SiStripFedCablingRcd'),
#                                                                    tag    = cms.string('SiStripFedCabling_GR10_v1_hlt')
#                                                                  )
#                                                        )
#                                     )
                                                                    
#process.load("CalibTracker.SiStripESProducers.SiStripConnectivity_cfi")


# substituting the specific payloads in teh GT
from CondCore.DBCommon.CondDBSetup_cfi import *
process.SetAPE = cms.ESSource("PoolDBESSource",CondDBSetup,
                              connect = cms.string('frontier://FrontierProd/CMS_COND_31X_STRIP'),
                              timetype = cms.string("runnumber"),
                              toGet = cms.VPSet(cms.PSet(record = cms.string('SiStripFedCablingRcd'),
                                                         tag = cms.string('SiStripFedCabling_GR10_v1_hlt')
                                                         )
                                                )
                              )
process.es_prefer_SetAPE = cms.ESPrefer("PoolDBESSource", "SetAPE")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string('/afs/cern.ch/cms/DB/conddb')
    ),
    timetype = cms.untracked.string('runnumber'),
    connect = cms.string('sqlite_file:dbfile.db'),
    toPut = cms.VPSet(cms.PSet(
        record = cms.string('SiStripApvGainRcd'),
        tag = cms.string('SiStripApvGainRcd_v1_hltvalidation')
    ), cms.PSet(
        record = cms.string('SiStripBadStripRcd'),
        tag = cms.string('SiStripBadChannel_v1_hltvalidation')
    ))
)

process.prod = cms.EDAnalyzer("SiStripApvGainFromFileBuilder",
    heightThreshold = cms.untracked.double(50.),
    doRecovery      = cms.untracked.bool(False),
    tickFile        = cms.untracked.FileInPath("CondTools/SiStrip/data/26May2015_timing.root"),
    referenceFile   = cms.untracked.FileInPath("CondTools/SiStrip/data/05May2015_timing.root"),#Oct2012_timing.root"),#11Apr2015_timing.root"),
    #recoveryList    = cms.untracked.string("CondTools/SiStrip/SiStripApvGainSummary.txt"), #CondTools/SiStrip/data/bad.txt"),
    outputMaps      = cms.untracked.bool(True),
    outputSummary   = cms.untracked.bool(True),
    putDummyIntoLowChannels = cms.untracked.bool(True),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('tickmark_analysis.root')
)

process.p = cms.Path(process.prod)


