#include "CondTools/SiStrip/plugins/SiStripLorentzAngleFromGainTree.h"

#include "CalibTracker/SiStripCommon/interface/SiStripDetInfoFileReader.h"
#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"

#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"


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
    SiStripGeometry geometry = (*el).second[0];
    if( geometry.isDummy() ) continue;
    uint32_t det_id = geometry.det_id; //(*el).first;
    //int structure_id = geometry.subdetectorStructure( det_id );

    //std::pair< int, std::string > structure = std::pair<int,std::string>(structure_id,geometry.subdetectorType(det_id) );
    int structure_id = geometry.subdetectorStructure();
    std::pair< int, std::string > structure = std::pair<int,std::string>(structure_id,geometry.subdetectorType() );

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
    //std::cout << "vector size at instertion: " << apvs.size() << ", APV number: " << nAPVs << std::endl;
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
                                                           edm::ESHandle<SiStripDetCabling>& cabling,
                                                           edm::ESHandle<TrackerGeometry>& tkGeometry)
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

      int subdetector,structure;
      geometry.whoIam(subdetector,structure);

      //global Orientation of local coordinate system of dets/detUnits
      const Surface& surface = tkGeometry->idToDet(geometry.det_id)->surface();
      const Surface::PositionType &gPModule = tkGeometry->idToDet(geometry.det_id)->position();

      if(subdetector==SiStripGeometry::TIB || subdetector==SiStripGeometry::TOB){
        LocalPoint  lWDirection(0.,0.,1.);
        GlobalPoint gWDirection = surface.toGlobal(lWDirection);
        float dRw = gWDirection.perp() - gPModule.perp();
        LocalPoint  lVDirection(0.,1.,0.);
        GlobalPoint gVDirection = surface.toGlobal(lVDirection);
        //float dRv = gVDirection.perp() - gPModule.perp();

        //LocalPoint  lUDirection(1.,0.,0.);
        //GlobalPoint gUDirection = surface.toGlobal(lUDirection);
        //float dRu = gUDirection.perp() - gPModule.perp();

        float dz = gVDirection.z() - gPModule.z();


        if ( dRw>0. && dz>0. ) geometry.localFrame = 1;
        if ( dRw>0. && dz<0. ) geometry.localFrame = 2;
        if ( dRw<0. && dz>0. ) geometry.localFrame = 3;
        if ( dRw<0. && dz<0. ) geometry.localFrame = 4;
/*
        if (geometry.localFrame==0) {
          std::cout << "dRw: " << dRw << ",  dRv: " << dRv << ",  dz: " << dz << std::endl;
          std::cout << "gVDirection.perp(): " << gVDirection.perp() << ",  gPModule.perp(): " << gPModule.perp() << std::endl;
          std::cout << "det_id: " << geometry.det_id
                    << ",  nAPVs: " << geometry.nAPVs
                    << ",  offlineAPV_id: " << geometry.offlineAPV_id
                    << ",  sensor: " << geometry.sensor << std::endl;
          std::cout << geometry.subdetectorType() << std::endl;

          std::cout << "FED_id: " << geometry.FED_id
                    << ",  FED_ch: " << geometry.FED_ch
                    << ",  i2cAdd: " << geometry.i2cAdd
                    << ",  FEC_crate: " << geometry.FEC_crate
                    << ",  FEC_slot: "  << geometry.FEC_slot
                    << ",  FEC_ring: " << geometry.FEC_ring
                    << ",  CCU_addr: " << geometry.CCU_addr
                    << ",  CCU_chan: " << geometry.CCU_chan << std::endl;
        }
*/
      }

      if(subdetector==SiStripGeometry::TID || subdetector==SiStripGeometry::TEC){
        //LocalPoint  lVDirection(0.,1.,0.);
        //GlobalPoint gVDirection = surface.toGlobal(lVDirection);
        //float dR = gVDirection.perp() - gPModule.perp();
        //if (dR>=0.) geometry.localDirection = 1;
        //else geometry.localDirection = -1;
        geometry.localFrame = 0;
      }



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
      //std::cout << "vector size at extraction: " << element.size() << ",  offline id: " << geometry.offlineAPV_id << std::endl;
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
  calibData_    = 0;

  moduleId_     = 0;
  localDirX_    = 0;
  localDirY_    = 0;
  localDirZ_    = 0;
  firstStrip_   = 0;
  clusterWidth_ = 0;
}


void
SiStripLorentzAngleFromGainTree::beginJob() {
  // Register to the TFileService
  edm::Service<TFileService> fs;
  TFileDirectory sensors = fs->mkdir("SensorType");

  // Book sensor type histograms
  for (unsigned int i=0; i<14; i++) {
    SiStripGeometry geometry;
    std::string tag         = "angle_vs_amplitude_"+geometry.sensorType( i+1 );
    std::string tit         = "Angle versus Amplitude for "+geometry.sensorType( i+1 );
    h_angle_vs_amplitude[i] = sensors.make<TH2F>(tag.c_str(),tit.c_str(),100,-1.5,1.5,10,0.5,11);
  }

  calibData_ = new TChain("gainCalibrationTree/tree");
  //calibData_ = new TChain("commonCalibrationTree/tree");

  //Collect the input calibration Trees
  edm::LogInfo("Workflow") << "@SUB=beginJob" << "Accessing input calibration file ..." << std::endl;
  std::vector< std::string >::const_iterator el = gainTreeList_.begin();
  while ( el!=gainTreeList_.end() ) {
    std::string filename = (*el);
    if ( filename.find("/store/")!=std::string::npos ) filename.insert(0,"root://cmsxrootd.fnal.gov//");

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
  iSetup.get<TrackerDigiGeometryRecord>().get(tkGeometry_);  


  siStripSetup_.initialize( reader, detCabling_, tkGeometry_);



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
      histo.substructure.push_back( (*el).first );
      std::string tag = (*el).second+"_ava";
      std::string tit = "Cluster amplitude versus track angle ";
      histo.h_angle_vs_amplitude.push_back( root_dir.make<TH2F>(tag.c_str(),tit.c_str(),100,-1.5,1.5,10,0.5,11) );
      ++el;
    }
    siStripHistos_.push_back( histo );
  }

  //Process the calibration Tree and fill histograms
  int Entries = calibData_->GetEntries();
  edm::LogInfo("Workflow") << "Processing " << Entries << " entries ..." << std::endl;

  TFile* current_file = 0;
  for(int i=0; i<=Entries;++i) {
    calibData_->GetEntry(i);

    if( current_file!=calibData_->GetCurrentFile() ) {

      current_file = calibData_->GetCurrentFile();

      if ( calibData_->GetBranch("GainCalibrationrawid")==0 ) throw("GainCalibrationrawid branch not found!");
      calibData_->SetBranchAddress("GainCalibrationrawid",&moduleId_);

      if ( calibData_->GetBranch("GainCalibrationlocaldirx")==0 ) throw("GainCalibrationlocaldirx");
      calibData_->SetBranchAddress("GainCalibrationlocaldirx",&localDirX_);

      if ( calibData_->GetBranch("GainCalibrationlocaldiry")==0 ) throw("GainCalibrationlocaldiry");
      calibData_->SetBranchAddress("GainCalibrationlocaldiry",&localDirY_);

      if ( calibData_->GetBranch("GainCalibrationlocaldirz")==0 ) throw("GainCalibrationlocaldirz");
      calibData_->SetBranchAddress("GainCalibrationlocaldirz",&localDirZ_);

      if ( calibData_->GetBranch("GainCalibrationfirststrip")==0 ) throw("GainCalibrationfirststrip");
      calibData_->SetBranchAddress("GainCalibrationfirststrip",&firstStrip_);

      if ( calibData_->GetBranch("GainCalibrationnstrips")==0 ) throw("GainCalibrationnstrips");
      calibData_->SetBranchAddress("GainCalibrationnstrips",&clusterWidth_);

      LogDebug("Debug") << "Opening " << calibData_->GetCurrentFile()->GetName() << std::endl;

    }

    for(unsigned int cluster=0;cluster<moduleId_->size();++cluster) {

      uint32_t det_id    = moduleId_->at(cluster);
      float directionX   = localDirX_->at(cluster);
      float directionY   = localDirY_->at(cluster);
      float directionZ   = localDirZ_->at(cluster);
      unsigned int strip = firstStrip_->at(cluster);
      unsigned int width = clusterWidth_->at(cluster);    

      float tanTheta =  directionX/directionZ;

      uint16_t offlineAPV_id = strip/128;
      //std::cout << "strip: " << strip << ",  offlineAPV_id: " << offlineAPV_id << std::endl;
      SiStripGeometry geometry = siStripSetup_.find( det_id, offlineAPV_id );

      int subdetector = 99;
      int structure   = 99;
      geometry.whoIam(subdetector,structure);
      //std::cout << subdetector << "  " << structure << std::endl;

      if( subdetector==SiStripGeometry::UNKNOWN ) continue;


      TH2F* h_per_sensor   = h_angle_vs_amplitude[geometry.sensor-1];
      TH2F* h_per_detector = 0;
      for(unsigned int substr=0;substr<siStripHistos_[subdetector].h_angle_vs_amplitude.size();++substr) {
        if( siStripHistos_[subdetector].substructure[substr]==structure) 
           h_per_detector = siStripHistos_[subdetector].h_angle_vs_amplitude[substr];
      }


      LogDebug("Debug") << "RawId: " << det_id
                        << "  clusterWidth: " << std::setw(2) << width
                        << "  dirX: " << std::setprecision(3) << std::setw(5) << directionX
                        << "  dirY: " << std::setprecision(3) << std::setw(5) << directionY 
                        << "  dirZ: " << std::setprecision(3) << std::setw(5) << directionZ 
                        << "  tan(theta): " << std::setprecision(3) << std::setw(5) << tanTheta << std::endl;

      h_per_sensor->Fill(tanTheta,width);
      h_per_detector->Fill(tanTheta,width);
    }
  }
}
