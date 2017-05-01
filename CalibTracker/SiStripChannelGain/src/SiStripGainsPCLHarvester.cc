// system include files
#include <memory>

// CMSSW includes
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
#include "CalibFormats/SiStripObjects/interface/SiStripGain.h"
#include "CalibFormats/SiStripObjects/interface/SiStripQuality.h"
#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
#include "CalibTracker/Records/interface/SiStripGainRcd.h"
#include "CalibTracker/Records/interface/SiStripQualityRcd.h"
#include "CondFormats/SiStripObjects/interface/SiStripApvGain.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"
#include "DataFormats/SiStripDetId/interface/SiStripSubStructure.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"

// user include files
#include "CalibTracker/SiStripChannelGain/interface/SiStripGainsPCLHarvester.h"
#include "CalibTracker/SiStripChannelGain/interface/APVGainHelpers.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include <iostream>
#include <sstream>

//********************************************************************************//
SiStripGainsPCLHarvester::SiStripGainsPCLHarvester(const edm::ParameterSet& ps):
  doStoreOnDB(false),
  GOOD(0),
  BAD(0),
  MASKED(0),
  NStripAPVs(0),
  NPixelDets(0),
  bareTkGeomPtr_(nullptr)
{

  m_Record                = ps.getUntrackedParameter<std::string> ("Record"             , "SiStripApvGainRcd");
  CalibrationLevel        = ps.getUntrackedParameter<int>         ("CalibrationLevel"   ,   0);
  MinNrEntries            = ps.getUntrackedParameter<double>      ("minNrEntries"       ,  20);
  m_DQMdir                = ps.getUntrackedParameter<std::string> ("DQMdir"             , "AlCaReco/SiStripGains");
  m_calibrationMode       = ps.getUntrackedParameter<std::string> ("calibrationMode"    , "StdBunch");
  tagCondition_NClusters  = ps.getUntrackedParameter<double>      ("NClustersForTagProd", 2E8);
  tagCondition_GoodFrac   = ps.getUntrackedParameter<double>      ("GoodFracForTagProd" , 0.95);
  doChargeMonitorPerPlane = ps.getUntrackedParameter<bool>        ("doChargeMonitorPerPlane" ,  false);
  VChargeHisto            = ps.getUntrackedParameter<std::vector<std::string> >  ("ChargeHisto");
  
  //Set the monitoring element tag and store
  dqm_tag_.reserve(7);
  dqm_tag_.clear();
  dqm_tag_.push_back( "StdBunch" );      // statistic collection from Standard Collision Bunch @ 3.8 T
  dqm_tag_.push_back( "StdBunch0T" );    // statistic collection from Standard Collision Bunch @ 0 T
  dqm_tag_.push_back( "AagBunch" );      // statistic collection from First Collision After Abort Gap @ 3.8 T
  dqm_tag_.push_back( "AagBunch0T" );    // statistic collection from First Collision After Abort Gap @ 0 T
  dqm_tag_.push_back( "IsoMuon" );       // statistic collection from Isolated Muon @ 3.8 T
  dqm_tag_.push_back( "IsoMuon0T" );     // statistic collection from Isolated Muon @ 0 T
  dqm_tag_.push_back( "Harvest" );       // statistic collection: Harvest

}

//********************************************************************************//
// ------------ method called for each event  ------------
void
SiStripGainsPCLHarvester::beginRun(edm::Run const& run, const edm::EventSetup& iSetup)
{
  using namespace edm;

  this->checkBookAPVColls(iSetup); // check whether APV colls are booked and do so if not yet done

  edm::ESHandle<SiStripGain> gainHandle;
  iSetup.get<SiStripGainRcd>().get(gainHandle);
  if(!gainHandle.isValid()){edm::LogError("SiStripGainPCLHarvester")<< "gainHandle is not valid\n"; exit(0);}

  edm::ESHandle<SiStripQuality> SiStripQuality_;
  iSetup.get<SiStripQualityRcd>().get(SiStripQuality_);

  for(unsigned int a=0;a<APVsCollOrdered.size();a++){

    std::shared_ptr<stAPVGain> APV = APVsCollOrdered[a];

    if(APV->SubDet==PixelSubdetector::PixelBarrel || APV->SubDet==PixelSubdetector::PixelEndcap) continue;
    
    APV->isMasked      = SiStripQuality_->IsApvBad(APV->DetId,APV->APVId);
    	  
    if(gainHandle->getNumberOfTags()!=2){edm::LogError("SiStripGainPCLHarvester")<< "NUMBER OF GAIN TAG IS EXPECTED TO BE 2\n";fflush(stdout);exit(0);};		   
    float newPreviousGain = gainHandle->getApvGain(APV->APVId,gainHandle->getRange(APV->DetId, 1),1);
    if(APV->PreviousGain!=1 and newPreviousGain!=APV->PreviousGain)edm::LogWarning("SiStripGainPCLHarvester")<< "WARNING: ParticleGain in the global tag changed\n";
    APV->PreviousGain = newPreviousGain;
    
    float newPreviousGainTick = gainHandle->getApvGain(APV->APVId,gainHandle->getRange(APV->DetId, 0),0);
    if(APV->PreviousGainTick!=1 and newPreviousGainTick!=APV->PreviousGainTick){
      edm::LogWarning("SiStripGainPCLHarvester")<< "WARNING: TickMarkGain in the global tag changed\n"<< std::endl
						 <<" APV->SubDet: "<< APV->SubDet << " APV->APVId:" << APV->APVId << std::endl
						 <<" APV->PreviousGainTick: "<<APV->PreviousGainTick<<" newPreviousGainTick: "<<newPreviousGainTick<<std::endl;
    }
    APV->PreviousGainTick = newPreviousGainTick;  	  
  }
}

//********************************************************************************//
void SiStripGainsPCLHarvester::dqmEndJob(DQMStore::IBooker& ibooker_, DQMStore::IGetter& igetter_) {

  edm::LogInfo("SiStripGainsPCLHarvester") << "Starting harvesting statistics" << std::endl;
    
  std::string DQM_dir = m_DQMdir;
   
  std::string stag =  *(std::find(dqm_tag_.begin(), dqm_tag_.end(),m_calibrationMode));
  if(stag.size()!=0 && stag[0]!='_') stag.insert(0,1,'_');

  std::string cvi      = DQM_dir + std::string("/Charge_Vs_Index") + stag;
     
  MonitorElement* Charge_Vs_Index           = igetter_.get(cvi.c_str());
  
  if (Charge_Vs_Index==0) {
    edm::LogError("SiStripGainsPCLHarvester") << "Harvesting: could not retrieve " << cvi.c_str()
					      << ", statistics will not be summed!" << std::endl;
  } else {
    edm::LogInfo("SiStripGainsPCLHarvester") << "Harvesting " 
					     << (Charge_Vs_Index)->getTH2S()->GetEntries() << " more clusters" << std::endl;
  }

  algoComputeMPVandGain(Charge_Vs_Index);
  std::unique_ptr<SiStripApvGain> theAPVGains = this->getNewObject(Charge_Vs_Index);

  // write out the APVGains record
  edm::Service<cond::service::PoolDBOutputService> poolDbService;
  
  if( doStoreOnDB ){
    if( poolDbService.isAvailable() )
      poolDbService->writeOne(theAPVGains.get(), poolDbService->currentTime(),m_Record);
    else
      throw std::runtime_error("PoolDBService required.");
  } else {
    edm::LogInfo("SiStripGainsPCLHarvester") << "Will not produce payload!" << std::endl;  
  }


  //Collect the statistics for monitoring and validation
  gainQualityMonitor(ibooker_,Charge_Vs_Index);

}

//********************************************************************************//
void
SiStripGainsPCLHarvester::gainQualityMonitor(DQMStore::IBooker& ibooker_, const MonitorElement* Charge_Vs_Index) const {
  ibooker_.setCurrentFolder("AlCaReco/SiStripGainsHarvesting/");


  std::vector<MonitorElement*> new_charge_histos;
  std::vector<std::pair<std::string,std::string>> cnames = APVGain::monHnames(VChargeHisto,doChargeMonitorPerPlane,"newG2");
  for (unsigned int i=0;i<cnames.size();i++){
    //std::string name = (cnames[i]).first + stag;
    new_charge_histos.push_back( ibooker_.book1DD( (cnames[i]).first, (cnames[i]).second.c_str(), 100   , 0. , 1000. ) );
  }


  MonitorElement* MPV_Vs_EtaTIB  = ibooker_.book2DD("MPV_Vs_EtaTIB" ,"MPV_Vs_EtaTIB" , 50, -3.0, 3.0, 300, 0, 600);
  MonitorElement* MPV_Vs_EtaTID  = ibooker_.book2DD("MPV_Vs_EtaTID" ,"MPV_Vs_EtaTID" , 50, -3.0, 3.0, 300, 0, 600);
  MonitorElement* MPV_Vs_EtaTOB  = ibooker_.book2DD("MPV_Vs_EtaTOB" ,"MPV_Vs_EtaTOB" , 50, -3.0, 3.0, 300, 0, 600);
  MonitorElement* MPV_Vs_EtaTEC  = ibooker_.book2DD("MPV_Vs_EtaTEC" ,"MPV_Vs_EtaTEC" , 50, -3.0, 3.0, 300, 0, 600);
  MonitorElement* MPV_Vs_EtaTEC1 = ibooker_.book2DD("MPV_Vs_EtaTEC1","MPV_Vs_EtaTEC1", 50, -3.0, 3.0, 300, 0, 600);
  MonitorElement* MPV_Vs_EtaTEC2 = ibooker_.book2DD("MPV_Vs_EtaTEC2","MPV_Vs_EtaTEC2", 50, -3.0, 3.0, 300, 0, 600);

  MonitorElement* MPV_Vs_PhiTIB  = ibooker_.book2DD("MPV_Vs_PhiTIB" ,"MPV_Vs_PhiTIB" , 50, -3.4, 3.4, 300, 0, 600);
  MonitorElement* MPV_Vs_PhiTID  = ibooker_.book2DD("MPV_Vs_PhiTID" ,"MPV_Vs_PhiTID" , 50, -3.4, 3.4, 300, 0, 600);
  MonitorElement* MPV_Vs_PhiTOB  = ibooker_.book2DD("MPV_Vs_PhiTOB" ,"MPV_Vs_PhiTOB" , 50, -3.4, 3.4, 300, 0, 600);
  MonitorElement* MPV_Vs_PhiTEC  = ibooker_.book2DD("MPV_Vs_PhiTEC" ,"MPV_Vs_PhiTEC" , 50, -3.4, 3.4, 300, 0, 600);
  MonitorElement* MPV_Vs_PhiTEC1 = ibooker_.book2DD("MPV_Vs_PhiTEC1","MPV_Vs_PhiTEC1", 50, -3.4, 3.4, 300, 0, 600);
  MonitorElement* MPV_Vs_PhiTEC2 = ibooker_.book2DD("MPV_Vs_PhiTEC2","MPV_Vs_PhiTEC2", 50, -3.4, 3.4, 300, 0, 600);

  MonitorElement* NoMPV          = ibooker_.book2DD("NoMPV"         ,"NoMPV"         ,350, -350, 350, 240, 0, 120);

  MonitorElement* Gains          = ibooker_.book1DD("Gains"         ,"Gains"         ,                300, 0, 2);
  MonitorElement* MPVs           = ibooker_.book1DD("MPVs"          ,"MPVs"          ,                300, 0, 600);
  MonitorElement* MPVs320        = ibooker_.book1DD("MPVs320"       ,"MPVs320"       ,                300, 0, 600);
  MonitorElement* MPVs500        = ibooker_.book1DD("MPVs500"       ,"MPVs500"       ,                300, 0, 600);
  MonitorElement* MPVsTIB        = ibooker_.book1DD("MPVsTIB"       ,"MPVsTIB"       ,                300, 0, 600);
  MonitorElement* MPVsTID        = ibooker_.book1DD("MPVsTID"       ,"MPVsTID"       ,                300, 0, 600);
  MonitorElement* MPVsTIDP       = ibooker_.book1DD("MPVsTIDP"      ,"MPVsTIDP"      ,                300, 0, 600);
  MonitorElement* MPVsTIDM       = ibooker_.book1DD("MPVsTIDM"      ,"MPVsTIDM"      ,                300, 0, 600);
  MonitorElement* MPVsTOB        = ibooker_.book1DD("MPVsTOB"       ,"MPVsTOB"       ,                300, 0, 600);
  MonitorElement* MPVsTEC        = ibooker_.book1DD("MPVsTEC"       ,"MPVsTEC"       ,                300, 0, 600);
  MonitorElement* MPVsTECP       = ibooker_.book1DD("MPVsTECP"      ,"MPVsTECP"      ,                300, 0, 600);
  MonitorElement* MPVsTECM       = ibooker_.book1DD("MPVsTECM"      ,"MPVsTECM"      ,                300, 0, 600);
  MonitorElement* MPVsTEC1       = ibooker_.book1DD("MPVsTEC1"      ,"MPVsTEC1"      ,                300, 0, 600);
  MonitorElement* MPVsTEC2       = ibooker_.book1DD("MPVsTEC2"      ,"MPVsTEC2"      ,                300, 0, 600);
  MonitorElement* MPVsTECP1      = ibooker_.book1DD("MPVsTECP1"     ,"MPVsTECP1"     ,                300, 0, 600);
  MonitorElement* MPVsTECP2      = ibooker_.book1DD("MPVsTECP2"     ,"MPVsTECP2"     ,                300, 0, 600);
  MonitorElement* MPVsTECM1      = ibooker_.book1DD("MPVsTECM1"     ,"MPVsTECM1"     ,                300, 0, 600);
  MonitorElement* MPVsTECM2      = ibooker_.book1DD("MPVsTECM2"     ,"MPVsTECM2"     ,                300, 0, 600);

  MonitorElement* MPVError       = ibooker_.book1DD("MPVError"      ,"MPVError"      ,                150, 0, 150);
  MonitorElement* MPVErrorVsMPV  = ibooker_.book2DD("MPVErrorVsMPV" ,"MPVErrorVsMPV" ,300,    0, 600, 150, 0, 150);
  MonitorElement* MPVErrorVsEta  = ibooker_.book2DD("MPVErrorVsEta" ,"MPVErrorVsEta" , 50, -3.0, 3.0, 150, 0, 150);
  MonitorElement* MPVErrorVsPhi  = ibooker_.book2DD("MPVErrorVsPhi" ,"MPVErrorVsPhi" , 50, -3.4, 3.4, 150, 0, 150);
  MonitorElement* MPVErrorVsN    = ibooker_.book2DD("MPVErrorVsN"   ,"MPVErrorVsN"   ,500,    0,1000, 150, 0, 150);

  MonitorElement* DiffWRTPrevGainTIB = ibooker_.book1DD("DiffWRTPrevGainTIB" ,"DiffWRTPrevGainTIB" , 250, 0,2);
  MonitorElement* DiffWRTPrevGainTID = ibooker_.book1DD("DiffWRTPrevGainTID" ,"DiffWRTPrevGainTID" , 250, 0,2);
  MonitorElement* DiffWRTPrevGainTOB = ibooker_.book1DD("DiffWRTPrevGainTOB" ,"DiffWRTPrevGainTOB" , 250, 0,2);
  MonitorElement* DiffWRTPrevGainTEC = ibooker_.book1DD("DiffWRTPrevGainTEC" ,"DiffWRTPrevGainTEC" , 250, 0,2);

  MonitorElement* GainVsPrevGainTIB  = ibooker_.book2DD("GainVsPrevGainTIB"  ,"GainVsPrevGainTIB"  , 100, 0,2, 100, 0,2);
  MonitorElement* GainVsPrevGainTID  = ibooker_.book2DD("GainVsPrevGainTID"  ,"GainVsPrevGainTID"  , 100, 0,2, 100, 0,2);
  MonitorElement* GainVsPrevGainTOB  = ibooker_.book2DD("GainVsPrevGainTOB"  ,"GainVsPrevGainTOB"  , 100, 0,2, 100, 0,2);
  MonitorElement* GainVsPrevGainTEC  = ibooker_.book2DD("GainVsPrevGainTEC"  ,"GainVsPrevGainTEC"  , 100, 0,2, 100, 0,2);


  for(unsigned int a=0;a<APVsCollOrdered.size();a++){

    std::shared_ptr<stAPVGain> APV = APVsCollOrdered[a];
    if(APV==NULL)continue;

    unsigned int  Index        = APV->Index;
    unsigned int  SubDet       = APV->SubDet;
    unsigned int  DetId        = APV->DetId;
    float         z            = APV->z;
    float         Eta          = APV->Eta;
    float         R            = APV->R;
    float         Phi          = APV->Phi;
    float         Thickness    = APV->Thickness;
    double        FitMPV       = APV->FitMPV;
    double        FitMPVErr    = APV->FitMPVErr;
    double        Gain         = APV->Gain;
    double        NEntries     = APV->NEntries;
    double        PreviousGain = APV->PreviousGain;


    std::vector<MonitorElement*> charge_histos = APVGain::FetchMonitor(new_charge_histos, DetId);
    TH2S *chvsidx = (Charge_Vs_Index)->getTH2S();
    int bin = chvsidx->GetXaxis()->FindBin(Index);
    TH1D* Proj = chvsidx->ProjectionY("proj",bin,bin);
    for (int binId=0; binId<Proj->GetXaxis()->GetNbins();binId++) {
      double new_charge = Proj->GetBinCenter(binId) / Gain;
      if (Proj->GetBinContent(binId)!=0.) {
        for (unsigned int h=0;h<charge_histos.size();h++) {
          TH1D* chisto = (charge_histos[h])->getTH1D();
          chisto->Fill(new_charge,Proj->GetBinContent(binId));
        }
      }
    }
    

    if(SubDet>=3 && FitMPV<0.) NoMPV->Fill(z,R);
    if(SubDet>=3 && FitMPV>=0.) {
      if(SubDet==3) MPV_Vs_EtaTIB->Fill(Eta,FitMPV);
      if(SubDet==4) MPV_Vs_EtaTID->Fill(Eta,FitMPV);
      if(SubDet==5) MPV_Vs_EtaTOB->Fill(Eta,FitMPV);
      if(SubDet==6) MPV_Vs_EtaTEC->Fill(Eta,FitMPV);
      if(SubDet==6 && Thickness<0.04) MPV_Vs_EtaTEC1->Fill(Eta,FitMPV);
      if(SubDet==6 && Thickness>0.04) MPV_Vs_EtaTEC2->Fill(Eta,FitMPV);
      if(SubDet==3) MPV_Vs_PhiTIB->Fill(Phi,FitMPV);
      if(SubDet==4) MPV_Vs_PhiTID->Fill(Phi,FitMPV);
      if(SubDet==5) MPV_Vs_PhiTOB->Fill(Phi,FitMPV);
      if(SubDet==6) MPV_Vs_PhiTEC->Fill(Phi,FitMPV);
      if(SubDet==6 && Thickness<0.04) MPV_Vs_PhiTEC1->Fill(Phi,FitMPV);
      if(SubDet==6 && Thickness>0.04) MPV_Vs_PhiTEC2->Fill(Phi,FitMPV);

      if(FitMPV>0.) Gains->Fill(Gain);
      MPVs->Fill(FitMPV);
      if(Thickness<0.04) MPVs320->Fill(Phi,FitMPV);
      if(Thickness>0.04) MPVs500->Fill(Phi,FitMPV);
      if(SubDet==3) MPVsTIB->Fill(FitMPV);
      if(SubDet==4) MPVsTID->Fill(FitMPV);
      if(SubDet==4 && Eta<0.) MPVsTIDM->Fill(FitMPV);
      if(SubDet==4 && Eta>0.) MPVsTIDP->Fill(FitMPV);
      if(SubDet==5) MPVsTOB->Fill(FitMPV);
      if(SubDet==6) MPVsTEC->Fill(FitMPV);
      if(SubDet==6 && Eta<0.) MPVsTECM->Fill(FitMPV);
      if(SubDet==6 && Eta>0.) MPVsTECP->Fill(FitMPV);
      if(SubDet==6 && Thickness<0.04) MPVsTEC1->Fill(FitMPV);
      if(SubDet==6 && Thickness>0.04) MPVsTEC2->Fill(FitMPV);
      if(SubDet==6 && Thickness<0.04 && Eta>0.) MPVsTECP1->Fill(FitMPV);
      if(SubDet==6 && Thickness>0.04 && Eta>0.) MPVsTECP2->Fill(FitMPV);
      if(SubDet==6 && Thickness<0.04 && Eta<0.) MPVsTECM1->Fill(FitMPV);
      if(SubDet==6 && Thickness>0.04 && Eta<0.) MPVsTECM2->Fill(FitMPV);

      MPVError->Fill(FitMPVErr);
      MPVErrorVsMPV->Fill(FitMPV,FitMPVErr);
      MPVErrorVsEta->Fill(Eta,FitMPVErr);
      MPVErrorVsPhi->Fill(Phi,FitMPVErr);
      MPVErrorVsN->Fill(NEntries,FitMPVErr);
    }

    if(SubDet==3 && PreviousGain!=0. ) DiffWRTPrevGainTIB->Fill(Gain/PreviousGain);
    if(SubDet==4 && PreviousGain!=0. ) DiffWRTPrevGainTID->Fill(Gain/PreviousGain);
    if(SubDet==5 && PreviousGain!=0. ) DiffWRTPrevGainTOB->Fill(Gain/PreviousGain);
    if(SubDet==6 && PreviousGain!=0. ) DiffWRTPrevGainTEC->Fill(Gain/PreviousGain);

    if(SubDet==3 ) GainVsPrevGainTIB->Fill(PreviousGain,Gain);
    if(SubDet==4 ) GainVsPrevGainTID->Fill(PreviousGain,Gain);
    if(SubDet==5 ) GainVsPrevGainTOB->Fill(PreviousGain,Gain);
    if(SubDet==6 ) GainVsPrevGainTEC->Fill(PreviousGain,Gain);

  }
}


//********************************************************************************//
void 
SiStripGainsPCLHarvester::algoComputeMPVandGain(const MonitorElement* Charge_Vs_Index) {

  unsigned int I=0;
  TH1F* Proj = NULL;
  double FitResults[6];
  double MPVmean = 300;

  if ( Charge_Vs_Index==0 ) {
    edm::LogError("SiStripGainsPCLHarvester") << "Harvesting: could not execute algoComputeMPVandGain method because "
					      << m_calibrationMode <<" statistics cannot be retrieved.\n"
					      << "Please check if input contains " 
					      << m_calibrationMode <<" data." << std::endl;
    return;
  }
  
  TH2S *chvsidx = (Charge_Vs_Index)->getTH2S();
 
  printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Fitting Charge Distribution  :");
  int TreeStep = APVsColl.size()/50;
  
  for(auto it = APVsColl.begin();it!=APVsColl.end();it++,I++){
    
    if(I%TreeStep==0){printf(".");fflush(stdout);}
    std::shared_ptr<stAPVGain> APV = it->second;
    if(APV->Bin<0) APV->Bin = chvsidx->GetXaxis()->FindBin(APV->Index);
    
    if(APV->isMasked){APV->Gain=APV->PreviousGain; MASKED++; continue;}
    
    Proj = (TH1F*)(chvsidx->ProjectionY("",chvsidx->GetXaxis()->FindBin(APV->Index),chvsidx->GetXaxis()->FindBin(APV->Index),"e"));
    if(!Proj)continue;
    
    if(CalibrationLevel==0){
    }else if(CalibrationLevel==1){
      int SecondAPVId = APV->APVId;
      if(SecondAPVId%2==0){    SecondAPVId = SecondAPVId+1; }else{ SecondAPVId = SecondAPVId-1; }
      std::shared_ptr<stAPVGain> APV2 = APVsColl[(APV->DetId<<4) | SecondAPVId];
      if(APV2->Bin<0) APV2->Bin = chvsidx->GetXaxis()->FindBin(APV2->Index);
      TH1F* Proj2 = (TH1F*)(chvsidx->ProjectionY("",APV2->Bin,APV2->Bin,"e"));
      if(Proj2){Proj->Add(Proj2,1);delete Proj2;}
    }else if(CalibrationLevel==2){
      for(unsigned int i=0;i<16;i++){  //loop up to 6APV for Strip and up to 16 for Pixels
	auto tmpit = APVsColl.find((APV->DetId<<4) | i); 
	if(tmpit==APVsColl.end())continue;
	std::shared_ptr<stAPVGain> APV2 = tmpit->second;
	if(APV2->DetId != APV->DetId || APV2->APVId == APV->APVId)continue;            
	if(APV2->Bin<0) APV2->Bin = chvsidx->GetXaxis()->FindBin(APV2->Index);
	TH1F* Proj2 = (TH1F*)(chvsidx->ProjectionY("",APV2->Bin,APV2->Bin,"e"));
	if(Proj2){Proj->Add(Proj2,1);delete Proj2;}
      }          
    } else {
      CalibrationLevel = 0;
      printf("Unknown Calibration Level, will assume %i\n",CalibrationLevel);
    }
    
    getPeakOfLandau(Proj,FitResults);
    APV->FitMPV      = FitResults[0];
    APV->FitMPVErr   = FitResults[1];
    APV->FitWidth    = FitResults[2];
    APV->FitWidthErr = FitResults[3];
    APV->FitChi2     = FitResults[4];
    APV->FitNorm     = FitResults[5];
    APV->NEntries    = Proj->GetEntries();
    
    if(IsGoodLandauFit(FitResults)){
      APV->Gain = APV->FitMPV / MPVmean;
      if(APV->SubDet>2)GOOD++;
    }else{
      APV->Gain = APV->PreviousGain;
      if(APV->SubDet>2)BAD++;
    }
    if(APV->Gain<=0)           APV->Gain  = 1;
    
    delete Proj;
  }printf("\n");
}

//********************************************************************************//
void 
SiStripGainsPCLHarvester::getPeakOfLandau(TH1* InputHisto, double* FitResults, double LowRange, double HighRange)
{ 
  FitResults[0]         = -0.5;  //MPV
  FitResults[1]         =  0;    //MPV error
  FitResults[2]         = -0.5;  //Width
  FitResults[3]         =  0;    //Width error
  FitResults[4]         = -0.5;  //Fit Chi2/NDF
  FitResults[5]         = 0;     //Normalization
  
  if( InputHisto->GetEntries() < MinNrEntries)return;
  
  // perform fit with standard landau
  TF1 MyLandau("MyLandau","landau",LowRange, HighRange);
  MyLandau.SetParameter(1,300);
  InputHisto->Fit(&MyLandau,"0QR WW");
  
  // MPV is parameter 1 (0=constant, 1=MPV, 2=Sigma)
  FitResults[0]         = MyLandau.GetParameter(1);  //MPV
  FitResults[1]         = MyLandau.GetParError(1);   //MPV error
  FitResults[2]         = MyLandau.GetParameter(2);  //Width
  FitResults[3]         = MyLandau.GetParError(2);   //Width error
  FitResults[4]         = MyLandau.GetChisquare() / MyLandau.GetNDF();  //Fit Chi2/NDF
  FitResults[5]         = MyLandau.GetParameter(0);
  
}

//********************************************************************************//
bool 
SiStripGainsPCLHarvester::IsGoodLandauFit(double* FitResults){
  if(FitResults[0] <= 0             )return false;
  //   if(FitResults[1] > MaxMPVError   )return false;
  //   if(FitResults[4] > MaxChi2OverNDF)return false;
  return true;   
}


//********************************************************************************//
// ------------ method called once each job just before starting event loop  ------------
void
SiStripGainsPCLHarvester::checkBookAPVColls(const edm::EventSetup& es){

  es.get<TrackerDigiGeometryRecord>().get( tkGeom_ );
  const TrackerGeometry *newBareTkGeomPtr = &(*tkGeom_);
  if (newBareTkGeomPtr == bareTkGeomPtr_) return; // already filled APVColls, nothing changed

  if (!bareTkGeomPtr_) { // pointer not yet set: called the first time => fill the APVColls
    auto const & Det = newBareTkGeomPtr->dets();
    
    unsigned int Index=0;

    for(unsigned int i=0;i<Det.size();i++){
      
      DetId  Detid  = Det[i]->geographicalId(); 
      int    SubDet = Detid.subdetId();
      
      if( SubDet == StripSubdetector::TIB ||  SubDet == StripSubdetector::TID ||
	  SubDet == StripSubdetector::TOB ||  SubDet == StripSubdetector::TEC  ){
	
	auto DetUnit     = dynamic_cast<const StripGeomDetUnit*> (Det[i]);
	if(!DetUnit)continue;
	
	const StripTopology& Topo     = DetUnit->specificTopology();	
	unsigned int         NAPV     = Topo.nstrips()/128;
	
	for(unsigned int j=0;j<NAPV;j++){
	  auto APV = std::make_shared<stAPVGain>();
	  APV->Index         = Index;
	  APV->Bin           = -1;
	  APV->DetId         = Detid.rawId();
	  APV->APVId         = j;
	  APV->SubDet        = SubDet;
	  APV->FitMPV        = -1;
	  APV->FitMPVErr     = -1;
	  APV->FitWidth      = -1;
	  APV->FitWidthErr   = -1;
	  APV->FitChi2       = -1;
	  APV->FitNorm       = -1;
	  APV->Gain          = -1;
	  APV->PreviousGain  = 1;
	  APV->PreviousGainTick  = 1;
	  APV->x             = DetUnit->position().basicVector().x();
	  APV->y             = DetUnit->position().basicVector().y();
	  APV->z             = DetUnit->position().basicVector().z();
	  APV->Eta           = DetUnit->position().basicVector().eta();
	  APV->Phi           = DetUnit->position().basicVector().phi();
	  APV->R             = DetUnit->position().basicVector().transverse();
	  APV->Thickness     = DetUnit->surface().bounds().thickness();
	  APV->NEntries	   = 0;
	  APV->isMasked      = false;
	  
	  APVsCollOrdered.push_back(APV);
	  APVsColl[(APV->DetId<<4) | APV->APVId] = APV;
	  Index++;
	  NStripAPVs++;
	} // loop on APVs
      } // if is Strips
    } // loop on dets
    
    for(unsigned int i=0;i<Det.size();i++){  //Make two loop such that the Pixel information is added at the end --> make transition simpler
      DetId  Detid  = Det[i]->geographicalId();
      int    SubDet = Detid.subdetId();
      if( SubDet == PixelSubdetector::PixelBarrel || PixelSubdetector::PixelEndcap ){
	auto DetUnit     = dynamic_cast<const PixelGeomDetUnit*> (Det[i]);
	if(!DetUnit) continue;
	
	const PixelTopology& Topo     = DetUnit->specificTopology();
	unsigned int         NROCRow  = Topo.nrows()/(80.);
	unsigned int         NROCCol  = Topo.ncolumns()/(52.);
	
	for(unsigned int j=0;j<NROCRow;j++){
	  for(unsigned int i=0;i<NROCCol;i++){
	    auto APV = std::make_shared<stAPVGain>();
	    APV->Index         = Index;
	    APV->Bin           = -1;
	    APV->DetId         = Detid.rawId();
	    APV->APVId         = (j<<3 | i);
	    APV->SubDet        = SubDet;
	    APV->FitMPV        = -1;
	    APV->FitMPVErr     = -1;
	    APV->FitWidth      = -1;
	    APV->FitWidthErr   = -1;
	    APV->FitChi2       = -1;
	    APV->Gain          = -1;
	    APV->PreviousGain  = 1;
	    APV->PreviousGainTick = 1;
	    APV->x             = DetUnit->position().basicVector().x();
	    APV->y             = DetUnit->position().basicVector().y();
	    APV->z             = DetUnit->position().basicVector().z();
	    APV->Eta           = DetUnit->position().basicVector().eta();
	    APV->Phi           = DetUnit->position().basicVector().phi();
	    APV->R             = DetUnit->position().basicVector().transverse();
	    APV->Thickness     = DetUnit->surface().bounds().thickness();
	    APV->isMasked      = false; //SiPixelQuality_->IsModuleBad(Detid.rawId());
	    APV->NEntries      = 0;
	    
	    APVsCollOrdered.push_back(APV);
	    APVsColl[(APV->DetId<<4) | APV->APVId] = APV;
	    Index++;
	    NPixelDets++;

	  } // loop on ROC cols
	} // loop on ROC rows
      } // if Pixel
    } // loop on Dets  
  }  //if (!bareTkGeomPtr_) ... 
  bareTkGeomPtr_ = newBareTkGeomPtr;
}

//********************************************************************************//
bool SiStripGainsPCLHarvester::produceTagFilter(const MonitorElement* Charge_Vs_Index){
  
  // The goal of this function is to check wether or not there is enough statistics
  // to produce a meaningful tag for the DB

  if( Charge_Vs_Index==0 ) {
    edm::LogError("SiStripGainsPCLHarvester") << "produceTagFilter -> Return false: could not retrieve the "
					      << m_calibrationMode <<" statistics.\n"
					      << "Please check if input contains " 
					      << m_calibrationMode << " data." << std::endl;
    return false; 
  }
  
  
  float integral = (Charge_Vs_Index)->getTH2S()->Integral();
  if( (Charge_Vs_Index)->getTH2S()->Integral(0,NStripAPVs+1, 0, 99999 ) < tagCondition_NClusters) {
    edm::LogWarning("SiStripGainsPCLHarvester") 
      << "calibrationMode  -> " << m_calibrationMode << "\n"
      << "produceTagFilter -> Return false: Statistics is too low : " << integral << std::endl;
    return false;
  }
  if((1.0 * GOOD) / (GOOD+BAD) < tagCondition_GoodFrac) {
    edm::LogWarning("SiStripGainsPCLHarvester")
      << "calibrationMode  -> " << m_calibrationMode << "\n"
      << "produceTagFilter ->  Return false: ratio of GOOD/TOTAL is too low: " << (1.0 * GOOD) / (GOOD+BAD) << std::endl;
    return false;
  }
  return true; 
}

//********************************************************************************//
std::unique_ptr<SiStripApvGain>
SiStripGainsPCLHarvester::getNewObject(const MonitorElement* Charge_Vs_Index) 
{
  std::unique_ptr<SiStripApvGain> obj = std::unique_ptr<SiStripApvGain>(new SiStripApvGain());
  
  if(!produceTagFilter(Charge_Vs_Index)){
    edm::LogWarning("SiStripGainsPCLHarvester")<< "getNewObject -> will not produce a paylaod because produceTagFilter returned false " << std::endl;       
    return obj;
  } else {
    doStoreOnDB=true;
  }
  
  std::vector<float> theSiStripVector;
  unsigned int PreviousDetId = 0; 
  for(unsigned int a=0;a<APVsCollOrdered.size();a++){
    std::shared_ptr<stAPVGain> APV = APVsCollOrdered[a];
    if(APV==NULL){ printf("Bug\n"); continue; }
    if(APV->SubDet<=2)continue;
    if(APV->DetId != PreviousDetId){
      if(!theSiStripVector.empty()){
	SiStripApvGain::Range range(theSiStripVector.begin(),theSiStripVector.end());
	if ( !obj->put(PreviousDetId,range) )  printf("Bug to put detId = %i\n",PreviousDetId);
      }
      theSiStripVector.clear();
      PreviousDetId = APV->DetId;
    }
    theSiStripVector.push_back(APV->Gain);
    
    LogDebug("SiStripGainsPCLHarvester")<<" DetId: "<<APV->DetId 
					<<" APV:   "<<APV->APVId
					<<" Gain:  "<<APV->Gain
					<<std::endl;

  }
  if(!theSiStripVector.empty()){
    SiStripApvGain::Range range(theSiStripVector.begin(),theSiStripVector.end());
    if ( !obj->put(PreviousDetId,range) )  printf("Bug to put detId = %i\n",PreviousDetId);
  }
  
  return obj;
}


//********************************************************************************//
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SiStripGainsPCLHarvester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
  
}

//********************************************************************************//
void SiStripGainsPCLHarvester::endRun(edm::Run const& run, edm::EventSetup const & isetup){
}

