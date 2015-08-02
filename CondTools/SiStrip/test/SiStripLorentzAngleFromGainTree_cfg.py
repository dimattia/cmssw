import FWCore.ParameterSet.Config as cms


process = cms.Process("ICALIB")
process.source = cms.Source("EmptyIOVSource",
    firstValue = cms.uint64(251252),#243664),#(243257), #(239640),
    lastValue = cms.uint64(251252),#243664), #(243257), #(239640),
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
process.GlobalTag = GlobalTag(process.GlobalTag, '74X_dataRun2_Express_v0', '') #''GR_P_V50', '')

#process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')


# substituting the specific payloads in teh GT
from CondCore.DBCommon.CondDBSetup_cfi import *
process.SetAPE = cms.ESSource("PoolDBESSource",CondDBSetup,
                              connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
                              timetype = cms.string("runnumber"),
                              toGet = cms.VPSet(cms.PSet(record = cms.string('SiStripFedCablingRcd'),
                                                         tag = cms.string('SiStripFedCabling_GR10_v1_hlt')
                                                         )
                                                )
                              )

process.es_prefer_SetAPE = cms.ESPrefer("PoolDBESSource", "SetAPE")


#process.myTrackerAlignmentErr = cms.ESSource("PoolDBESSource",CondDBSetup,
#                              connect = cms.string('frontier://FrontierProd/CMS_COND_31X_FROM21X'),
#                              timetype = cms.string("runnumber"),
#                              toGet = cms.VPSet(cms.PSet(record = cms.string('TrackerAlignmentErrorExtendedRcd'),
#                                                         tag = cms.string('TrackerIdealGeometryErrors210_mc')
#                                                        )
#                                               )
#                                            )

#process.es_prefer_myTrackerAlignmentErr = cms.ESPrefer("PoolDBESSource", "myTrackerAlignmentErr")


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
        #"/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR15/calibTree_251252.root",
        #"/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR15/calibTree_251252_100.root",
        #"/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR15/calibTree_251252_150.root",
        #"/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR15/calibTree_251252_200.root",
        #"/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR15/calibTree_251252_250.root",
        #"/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR15/calibTree_251252_350.root",
        #"/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR15/calibTree_251252_400.root",
        #"/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR15/calibTree_251252_50.root",
                                      
#        "/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR15/calibTree_251024.root",
#        "/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR15/calibTree_251025.root",
       # "/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR15/calibTree_251027.root",
       # "/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR12/calibTree_207905.root"
#         "calibTree_207905.root"
         "/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR15/calibTree_250866.root",
         "/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR15/calibTree_250867.root",
         "/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR15/calibTree_250868.root",
         "/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR15/calibTree_250869.root"
                                            )
)

process.TFileService = cms.Service("TFileService",
    #fileName = cms.string('lorentz_angle_new_run207295.root')#51252.root')
    fileName = cms.string('peak_0T_2015.root')
)

process.p = cms.Path(process.prod)


