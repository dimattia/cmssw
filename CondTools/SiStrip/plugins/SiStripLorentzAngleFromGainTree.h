#ifndef SiStripLorentzAngleFromGainTree_H
#define SiStripLorentzAngleFromGainTree_H

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"


#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#include <map>
#include <vector>
#include <sstream>
#include <stdint.h>


class TrackerGeometry;
class SiStripDetCabling;
class SiStripDetInfoFileReader;
class TH2F;
class TChain;

class SiStripLorentzAngleFromGainTree : public edm::EDAnalyzer {

 public:

  /** Brief: Stores geometry and cabling data per detector id.
   */
  struct SiStripGeometry {
                   typedef enum { TIB, TID, TOB, TEC, ALL, UNKNOWN } Subdetectors;
                   uint32_t  det_id;
                   uint16_t  nAPVs;
                   uint16_t  offlineAPV_id;
                   uint16_t  sensor;

                   int       localFrame;
                   
                   int       FED_id;
                   int       FED_ch;
                   int       i2cAdd;
                   int       FEC_crate;
                   int       FEC_slot;
                   int       FEC_ring;
                   int       CCU_addr;
                   int       CCU_chan;

                   bool      is_connected;

                   float     strip_length;

                   SiStripGeometry() {
                     det_id = 0xffffffff;
                     nAPVs  = 0xffff;
                     sensor = 0xffff;

                     localFrame = 0;

                     FED_id    = -1;
                     FED_ch    = -1;
                     i2cAdd    = -1;
                     FEC_crate = -1;
                     FEC_slot  = -1;
                     FEC_ring  = -1;
                     CCU_addr  = -1;
                     CCU_chan  = -1;

                     is_connected = false;

                     strip_length = 0.;
                   }

                   /** Brief Checks if SiStripGeometry is filled.
                    */
                   bool isDummy() { return det_id==0xffffffff; }  

                   /** Brief Returns the detector type.
                    */
                   void whoIam( int& subdetector , int& structure ) { 
                     if( isTIB(det_id) ) subdetector = TIB;
                     else if( isTID(det_id) ) subdetector = TID;
                     else if( isTOB(det_id) ) subdetector = TOB;
                     else if( isTEC(det_id) ) subdetector = TEC;  
                     else subdetector = UNKNOWN;

                     structure = subdetectorStructure( );
                   }

                   /** Brief: Returns the name of the sensor type
                    */
                   std::string sensorType(uint16_t Sensor) {
                     std::string types[15] = {"UNKNOWN",
                                              "IB1","IB2","OB1","OB2","W1A","W2A","W3A","W1B","W2B","W3B","W4","W5","W6","W7"};
                     if ( Sensor>15 ) return types[0];
                     return types[Sensor];
                   }

                   /** Brief Return the APV online Id.
                    */
                   uint16_t onlineAPV_id(void) { return i2cAdd%32; } 

                   /** Brief Return the APV offline Id.
                    */
                   uint16_t online2offline(uint16_t onlineAPV_id) 
                     { return (onlineAPV_id>=nAPVs)? onlineAPV_id -2 : onlineAPV_id; } 

                   /** Brief Check is det_id belogs to TIB.
                    */
                   bool isTIB(uint32_t DetId) { return ( ((DetId>>25)&0x7)==3 ); }

                   /** Brief Check is det_id belogs to TID.
                    */
                   bool isTID(uint32_t DetId) { return ( ((DetId>>25)&0x7)==4 );  }

                   /** Brief Check is det_id belogs to TOB.
                    */
                   bool isTOB(uint32_t DetId) { return ( ((DetId>>25)&0x7)==5 );  }

                   /** Brief Check is det_id belogs to TEC.
                    */
                   bool isTEC(uint32_t DetId) { return ( ((DetId>>25)&0x7)==6 );  }

                   /** Brief Return the subdetector structure (Layer or Wheel).
                    */
                   int subdetectorStructure() { 
                     if (isDummy()) return 0;

                     int st = 0;
                     if( isTIB(det_id) )      st = TIBDetId( det_id ).layer(); 
                     else if( isTID(det_id) ) st = -1*(2*(int)TIDDetId( det_id ).side()-3)*(int)TIDDetId( det_id ).wheel();
                     else if( isTOB(det_id) ) st = TOBDetId( det_id ).layer();
                     else if( isTEC(det_id) ) st = -1*(2*(int)TECDetId( det_id ).side()-3)*(int)TECDetId( det_id ).wheel();

                     int sgn = (st<0.)? -1 : 1;

                     return  st + sgn*localFrame*100;
                   }

                   /** Brief Check if the subdetector module is in positive side.
                    */ 
                   bool isZPlusSide(uint32_t DetId) {
                     if( isTIB(DetId) )      return TIBDetId( DetId ).isZPlusSide();
                     else if( isTID(DetId) ) return TIDDetId( DetId ).isZPlusSide();
                     else if( isTOB(DetId) ) return TOBDetId( DetId ).isZPlusSide();
                     else if( isTEC(DetId) ) return TECDetId( DetId ).isZPlusSide();
                     return 0;
                   }

                   /** Brief Check if the subdetector module is in positive side.
                    */
                   bool isZMinusSide(uint32_t DetId) {
                     if( isTIB(DetId) )      return TIBDetId( DetId ).isZMinusSide();
                     else if( isTID(DetId) ) return TIDDetId( DetId ).isZMinusSide();
                     else if( isTOB(DetId) ) return TOBDetId( DetId ).isZMinusSide();
                     else if( isTEC(DetId) ) return TECDetId( DetId ).isZMinusSide();
                     return 0;
                   }


                   /** Brief Return a string tag
                    */
                   std::string subdetectorType() {
                     if(isDummy()) return "";

                     int subDetId = ((det_id>>25)&0x7);
                     int subStrId = ( (int)std::fabs( subdetectorStructure() )) % 100;
                     bool negativeSide = isZMinusSide(det_id);
                     std::stringstream out;

                     switch(subDetId){
                     case 3: //TIB
                       out << "TIB_layer" << subStrId << "_noSide";
                       break;
                     case 4: //TID
                       out << "TID_wheel" << subStrId;
                       out << ( (negativeSide)? "_negativeSide" : "_positiveSide" );
                       break;
                     case 5: //TOB
                       out << "TOB_layer" << subStrId << "_noSide";
                       break;
                     case 6: //TEC
                       out << "TEC_wheel" << subStrId;
                       out << ( (negativeSide)? "_negativeSide" : "_positiveSide" );
                       break;
                     }

                     out << "_frame" << localFrame ;
                     return out.str();
                   }  

                 };

  /** Brief Manages the detector geometry and cabling setup.
   */
  class DetectorSetup {
    public:
    typedef std::map<uint32_t, std::vector<SiStripGeometry>> GeoMap;

    DetectorSetup();
    ~DetectorSetup();

    /** Brief Gathers the geometry for a given det_id and APV id.
     */
    SiStripGeometry find(uint32_t det_id, uint16_t offlineAPV_id=0) const;

    /** Brief Returns the subtructures for a given subdetector.
     */
    std::set< std::pair<int,std::string> > substructures( SiStripGeometry::Subdetectors det ) const;


    /** Brief Initialize the geometry and cabling map.
     */  
    void initialize(SiStripDetInfoFileReader& reader, 
                    edm::ESHandle<SiStripDetCabling>& cabling,
                    edm::ESHandle<TrackerGeometry>& tkGeometry);


    private:
    GeoMap detectorSetup_;  /*!< store the detector geometry and cabling setup. */

    /** Brief Return the element mapped by the detector id.
     */ 
    std::vector<SiStripGeometry>& retrieve_element(uint32_t det_id, unsigned int nAPVs); 
  };  




  /** Brief Constructor.
   */ 
  explicit SiStripLorentzAngleFromGainTree( const edm::ParameterSet& iConfig);

  /** Brief Destructor performing the memory cleanup.
   */ 
  ~SiStripLorentzAngleFromGainTree();

  /** Brief One dummy-event analysis to create the Lorentz Angle histograms.
   * This method creates the geometry structures, then loop over the gain trees
   * to fill the Lorentz Angle histograms. 
   */ 
  virtual void analyze(const edm::Event& , const edm::EventSetup& );

  /** Brief Routine to be called when the Job start.
   */ 
  virtual void beginJob();

  /** Brief Routine to be called when the Job end.
   */  
  virtual void endJob();


 private:
  struct DetectorHistograms {
           std::vector<int> substructure;
           std::vector<TH2F*> h_angle_vs_amplitude; /*!< cluster amplitude versus track angle. */
         };

  edm::FileInPath gfp_;                             /*!< File Path for the ideal geometry. */
  std::vector<std::string> gainTreeList_;           /*!< List of calibration tree to be processed. */

  TChain* calibData_;                               /*!< Contains the chained input calibration ntuple. */

  TH2F* h_angle_vs_amplitude[14];                   /*!< amplitude versus track angle per sensor type. */
  std::vector<DetectorHistograms> siStripHistos_;   /*!< histogram per detector type. */

  edm::ESHandle<SiStripDetCabling> detCabling_;     /*!< Handles the SiStrip cabling retrieved from the eventSetup. */
  edm::ESHandle<TrackerGeometry> tkGeometry_;       /*!< Handles the tracker geometry retrieved fromm the eventSetup. */
  DetectorSetup siStripSetup_;                      /*!< Describes the SiStrip geometry and cabling setup. */

  std::vector<uint32_t>* moduleId_;                 /*!< Detector Id from gain calibration ntuple. */
  std::vector<float>* localDirX_;                   /*!< Track direction along X from gain calibration ntuple. */
  std::vector<float>* localDirY_;                   /*!< Track direction along Y from gain calibration ntuple. */
  std::vector<float>* localDirZ_;                   /*!< Track direction along Z from gain calibration ntuple. */
  std::vector<unsigned short>* firstStrip_;         /*!< Index of the first strip in the cluster. */
  std::vector<unsigned short>* clusterWidth_;       /*!< Cluster amplitude from gain calibration ntuple. */
};

#endif
