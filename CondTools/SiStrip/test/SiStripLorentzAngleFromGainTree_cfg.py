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
#process.PoolDBOutputService = cms.Service("PoolDBOutputService",
#    BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
#    DBParameters = cms.PSet(
#        authenticationPath = cms.untracked.string('/afs/cern.ch/cms/DB/conddb')
#    ),
#    timetype = cms.untracked.string('runnumber'),
#    connect = cms.string('sqlite_file:dbfile.db'),
#    toPut = cms.VPSet(cms.PSet(
#        record = cms.string('SiStripApvGainRcd'),
#        tag = cms.string('SiStripApvGainRcd_v1_hltvalidation')
#    ), cms.PSet(
#        record = cms.string('SiStripBadStripRcd'),
#        tag = cms.string('SiStripBadChannel_v1_hltvalidation')
#    ))
#)

process.prod = cms.EDAnalyzer("SiStripLorentzAngleFromGainTree",
    inputCalibNtuple = cms.untracked.vstring(
        "/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR15/calibTree_251024.root",
        "/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR15/calibTree_251025.root",
        "/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR15/calibTree_251027.root",
                                            )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('lorentz_angle.root')
)

process.p = cms.Path(process.prod)


