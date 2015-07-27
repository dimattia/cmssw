#include "CondTools/SiStrip/plugins/SiStripLorentzAngleFromGainTree.h"
#include "CalibTracker/SiStripCommon/interface/SiStripDetInfoFileReader.h"
#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"

#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TChain.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <utility>
#include <algorithm>

#include <stdexcept>

#include <sys/types.h>
#include <sys/stat.h>


SiStripLorentzAngleFromGainTree::DetectorSetup::DetectorSetup() { }
SiStripLorentzAngleFromGainTree::DetectorSetup::~DetectorSetup() { }

std::set< std::pair<int, std::string> >
SiStripLorentzAngleFromGainTree::DetectorSetup::substructures( SiStripGeometry::Subdetectors detector ) const
{

  std::set< std::pair<int,std::string> > structureIds;

  for( GeoMap::const_iterator el= detectorSetup_.begin(); el!=detectorSetup_.end(); ++el ) {
    SiStripGeometry geometry;
    uint32_t det_id = (*el).first;
    int structure_id = geometry.subdetectorStructure( det_id );

    std::pair< int, std::string > structure = std::pair<int,std::string>(structure_id,geometry.subdetectorType(det_id) );

    if ( detector==SiStripGeometry::TIB && geometry.isTIB( det_id ) ) structureIds.insert( structure );
    else if ( detector==SiStripGeometry::TID && geometry.isTID( det_id ) ) structureIds.insert( structure );
    else if ( detector==SiStripGeometry::TOB && geometry.isTOB( det_id ) ) structureIds.insert( structure );
    else if ( detector==SiStripGeometry::TEC && geometry.isTEC( det_id ) ) structureIds.insert( structure );
  }
  return structureIds;
}

std::vector<SiStripLorentzAngleFromGainTree::SiStripGeometry>& 
SiStripLorentzAngleFromGainTree::DetectorSetup::retrieve_element(uint32_t det_id, unsigned int nAPVs)
{
  GeoMap::iterator element = detectorSetup_.find( det_id );

  // id det id not found, puts a dummy one then returns the associated empty vector
  if ( element==detectorSetup_.end() ) {
    std::vector<SiStripGeometry> apvs(nAPVs,SiStripGeometry());
    std::cout << "vector size at instertion: " << apvs.size() << ", APV number: " << nAPVs << std::endl;
    std::pair<GeoMap::iterator,bool> ret=detectorSetup_.insert(std::pair<uint32_t,std::vector<SiStripGeometry>>(det_id, apvs));
    if ( ret.second == false ) {
      edm::LogError("MapError") << "@SUB=DetectorSetup::retrieve_element" << "Cannot not insert geometry for detector id "
                                << det_id << " into the internal map: detector id already in the map." << std::endl;
      throw("Inconsistency when filling the geometry map.");
    }
    return retrieve_element(det_id,nAPVs);
  }

  return (*element).second;
}

void
SiStripLorentzAngleFromGainTree::DetectorSetup::initialize(SiStripDetInfoFileReader& reader,
                                                           edm::ESHandle<SiStripDetCabling>& cabling)
{
  detectorSetup_.clear();

  const std::map<uint32_t, SiStripDetInfoFileReader::DetInfo > DetInfos  = reader.getAllData();

  for(std::map<uint32_t, SiStripDetInfoFileReader::DetInfo >::const_iterator it = DetInfos.begin(); it != DetInfos.end(); it++){  
    // check if det id is correct and if it is actually cabled in the detector
    if( it->first==0 || it->first==0xFFFFFFFF ) {
      edm::LogError("DetIdNotGood") << "@SUB=DetectorSetup::initialize" << "Wrong det id: " << it->first
                                    << "  ... neglecting!" << std::endl;
      continue;
    }


    // For the cabled det_id retrieve the number of APV connected
    // to the module with the FED cabling
    uint16_t nAPVs = 0;
    const std::vector<const FedChannelConnection*> connection = cabling->getConnections(it->first);

    for (unsigned int ca = 0; ca<connection.size(); ca++) {
      if ( connection[ca]!=0 )  {
        nAPVs+=( connection[ca] )->nApvs();
        break;
      }
    }

    // check consistency among FED cabling and ideal cabling, exit on error
    if (connection.size()!=0 && nAPVs != (uint16_t)it->second.nApvs) {
      edm::LogError("SiStripCablingError") << "@SUB=DetectorSetup::initialize"
           << "det id " << it->first << ": APV number from FedCabling (" << nAPVs
           << ") is different from the APV number retrieved from the ideal cabling ("
           << it->second.nApvs << ")." << std::endl;
      throw("Inconsistency on the number of APVs.");
    }

    // Fill the gain in the DB object, apply the logic for the dummy values
    for(unsigned short j=0; j<it->second.nApvs; j++){
      SiStripGeometry geometry;
      geometry.det_id = it->first;
      geometry.offlineAPV_id = j;
      geometry.nAPVs  = it->second.nApvs;
      geometry.sensor = SiStripDetId(it->first).moduleGeometry();

      // put dummy values for the cabling; to be overwritten if module is actually cabled
      geometry.FED_id       = -1;
      geometry.FED_ch       = -1;
      geometry.i2cAdd       = -1;
      geometry.FEC_crate    = -1;
      geometry.FEC_slot     = -1;
      geometry.FEC_ring     = -1;
      geometry.CCU_addr     = -1;
      geometry.CCU_chan     = -1;
      geometry.is_connected = false;

      geometry.strip_length = reader.getNumberOfApvsAndStripLength(it->first).second;

      // fill the cabling for connected modules
      for (unsigned int ca = 0; ca<connection.size(); ca++) {
        if( connection[ca]!=0 && ( geometry.online2offline( (connection[ca])->i2cAddr(j%2)%32 )==geometry.offlineAPV_id )) {
          geometry.is_connected = (connection[ca])->isConnected();
          geometry.FED_id = (connection[ca])->fedId();
          geometry.FED_ch = (connection[ca])->fedCh();
          geometry.i2cAdd = (connection[ca])->i2cAddr( j%2 );
          geometry.FEC_crate = (connection[ca])->fecCrate();
          geometry.FEC_slot  = (connection[ca])->fecSlot();
          geometry.FEC_ring  = (connection[ca])->fecRing();
          geometry.CCU_addr  = (connection[ca])->ccuAddr();
          geometry.CCU_chan  = (connection[ca])->ccuChan();
        }
      }

      // Fill the geometry setup map
      std::vector<SiStripGeometry>& element = retrieve_element( it->first, it->second.nApvs );
      std::cout << "vector size at extraction: " << element.size() << ",  offline id: " << geometry.offlineAPV_id << std::endl;
      element[ geometry.offlineAPV_id ] = geometry;

    } // end loop on APVs

  } // end loop on det ids
}

SiStripLorentzAngleFromGainTree::SiStripGeometry
SiStripLorentzAngleFromGainTree::DetectorSetup::find( uint32_t det_id, uint16_t offlineAPV_id) const {
  SiStripGeometry geometry = SiStripGeometry();
  GeoMap::const_iterator element;
  element = detectorSetup_.find( det_id );
  if ( element!=detectorSetup_.end() ) geometry = (*element).second.at(offlineAPV_id);
  return geometry;
}

SiStripLorentzAngleFromGainTree::~SiStripLorentzAngleFromGainTree() { if(calibData_) delete calibData_; }

SiStripLorentzAngleFromGainTree::SiStripLorentzAngleFromGainTree( const edm::ParameterSet& iConfig ):
  gfp_(iConfig.getUntrackedParameter<edm::FileInPath>("geoFile",edm::FileInPath("CalibTracker/SiStripCommon/data/SiStripDetInfo.dat"))),
  gainTreeList_(iConfig.getUntrackedParameter<std::vector<std::string>>("inputCalibNtuple"))
{
  calibData_ = 0;
}


void
SiStripLorentzAngleFromGainTree::beginJob() {
  // Register to the TFileService
  edm::Service<TFileService> fs;
  TFileDirectory sensors = fs->mkdir("SensorType");

  // Book sensor type histograms
  for (unsigned int i=0; i<14; i++) {
    SiStripGeometry geometry;
    std::string tag         = "angle_vs_amplitude_"+geometry.sensorType( i );
    std::string tit         = "Angle versus Amplitude for "+geometry.sensorType( i );
    h_angle_vs_amplitude[i] = sensors.make<TH2F>(tag.c_str(),tit.c_str(),100,-2.,2.,20,0.5,21);
  }

  calibData_ = new TChain("gainCalibrationTree/tree");

  //Collect the input calibration Trees
  edm::LogInfo("Workflow") << "@SUB=beginJob" << "Accessing input calibration file ..." << std::endl;
  std::vector< std::string >::const_iterator el = gainTreeList_.begin();
  while ( el!=gainTreeList_.end() ) {
    std::string filename = (*el);
    if ( filename.find("cmsxrootd")==std::string::npos ) filename.insert(0,"root://cmsxrootd.fnal.gov//");

    edm::LogVerbatim("Workflow") << " ... " << filename << std::endl;

    calibData_->Add( filename.c_str() );
    
    ++el;
  }
}

void
SiStripLorentzAngleFromGainTree::endJob() {
}


void SiStripLorentzAngleFromGainTree::analyze(const edm::Event& evt, const edm::EventSetup& iSetup){

  // Retrieve the SiStripDetCabling description
  iSetup.get<SiStripDetCablingRcd>().get( detCabling_ );

  // Gather the geometry description
  SiStripDetInfoFileReader reader(gfp_.fullPath());

  siStripSetup_.initialize( reader, detCabling_ );


  // instanciate the histograms per detector structure
  edm::Service<TFileService> fs;
  std::string dirname[4] = {"TIB","TID","TOB","TEC"};
  for (int subDetector=SiStripGeometry::TIB;subDetector<=SiStripGeometry::TEC;++subDetector) {

    TFileDirectory root_dir = fs->mkdir( dirname[subDetector] );
    DetectorHistograms histo;

    SiStripGeometry::Subdetectors s = static_cast<SiStripGeometry::Subdetectors>(subDetector);
    std::set< std::pair<int,std::string> > subStructure = siStripSetup_.substructures( s );
    std::set< std::pair<int,std::string> >::const_iterator el = subStructure.begin();

    while ( el!=subStructure.end() ) {
      histo.substructure = (*el).first;
      std::string tag = (*el).second+"_ava";
      std::string tit = "Cluster amplitude versus track angle ";
      histo.h_angle_vs_amplitude.push_back( root_dir.make<TH2F>(tag.c_str(),tit.c_str(),100,-2.,2.,20,0.5,21) );
      ++el;
    }
  }

  //Process the calibration Tree and fill histograms
  int Entries = calibData_->GetEntries();
  LogDebug("Debug") << "Entries in calibration ntuple: " << Entries << std::endl;
  for(int i=0; i<=Entries;++i) {
  }
}
