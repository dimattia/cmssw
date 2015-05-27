#include "CondTools/SiStrip/plugins/SiStripApvGainFromFileBuilder.h"
#include "CalibTracker/SiStripCommon/interface/SiStripDetInfoFileReader.h"
#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"

#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <utility>
#include <algorithm>

#include <stdexcept>

#include <sys/types.h>
#include <sys/stat.h>




struct clean_up {
   void operator()(SiStripApvGainFromFileBuilder::Gain* el) {
       if (el!=0) {
           el->clear();
           delete el;
           el = 0;
       }
   }
} CleanUp;


bool ends_with(std::string const & value, std::string const & ending) {
    if( ending.size() > value.size() ) return false;
    return std::equal( ending.rbegin(), ending.rend(), value.rbegin() );
}


SiStripApvGainFromFileBuilder::~SiStripApvGainFromFileBuilder() {
    for_each(gains_.begin(), gains_.end(), CleanUp);
    for_each(negative_gains_.begin(), negative_gains_.end(), CleanUp);
    for_each(null_gains_.begin(), null_gains_.end(), CleanUp);

    if (siStripQuality_!=0) delete siStripQuality_;
}

SiStripApvGainFromFileBuilder::SiStripApvGainFromFileBuilder( const edm::ParameterSet& iConfig ):
  gfp_(iConfig.getUntrackedParameter<edm::FileInPath>("geoFile",edm::FileInPath("CalibTracker/SiStripCommon/data/SiStripDetInfo.dat"))),
  tfp_(iConfig.getUntrackedParameter<edm::FileInPath>("tickFile",edm::FileInPath("CondTools/SiStrip/data/tickheight.txt"))),
  rfp_(iConfig.getUntrackedParameter<edm::FileInPath>("referenceFile",edm::FileInPath("CondTools/SiStrip/data/tickheight.txt"))),
  nfp_(iConfig.getUntrackedParameter<edm::FileInPath>("noiseFile",edm::FileInPath("CondTools/SiStrip/data/noise_May.root"))),
  recoveryList_(iConfig.getUntrackedParameter<edm::FileInPath>("recoveryList",edm::FileInPath("CondTools/SiStrip/data/bad.txt"))),
  heightThreshold_(iConfig.getUntrackedParameter<double>("heightThreshold",0.)), 
  goodHeightLimit_(iConfig.getUntrackedParameter<double>("goodHeightLimit",440.)),
  badHeightLimit_(iConfig.getUntrackedParameter<double>("badHeightLimit",400.)),
  dummyAPVGain_(iConfig.getUntrackedParameter<double>("dummyAPVGain",690./640.)),
  doRecovery_(iConfig.getUntrackedParameter<bool>("doRecovery",1)),
  putDummyIntoUncabled_(iConfig.getUntrackedParameter<bool>("putDummyIntoUncabled",1)),
  putDummyIntoUnscanned_(iConfig.getUntrackedParameter<bool>("putDummyIntoUnscanned",1)),
  putDummyIntoOffChannels_(iConfig.getUntrackedParameter<bool>("putDummyIntoOffChannels",1.)),
  putDummyIntoBadChannels_(iConfig.getUntrackedParameter<bool>("putDummyIntoBadChannels",1.)),
  putDummyIntoLowChannels_(iConfig.getUntrackedParameter<bool>("putDummyIntoLowChannels",0.)),
  outputMaps_(iConfig.getUntrackedParameter<bool>("outputMaps",0)),
  outputSummary_(iConfig.getUntrackedParameter<bool>("outputSummary",0))
{
  siStripQuality_ = 0;
  siStripQuality_ = new SiStripQuality();
}


void
SiStripApvGainFromFileBuilder::beginJob() {
  // Register to the TFileService
  edm::Service<TFileService> fs;

  // Book histogram
  h_tickmark_height_0 = fs->make<TH1F>("tickmark_height_0", "Tickmark height with AOH gain at 0", 110,0.,1100.);
  h_tickmark_height_1 = fs->make<TH1F>("tickmark_height_1", "Tickmark height with AOH gain at 1", 110,0.,1100.);
  h_tickmark_height_2 = fs->make<TH1F>("tickmark_height_2", "Tickmark height with AOH gain at 2", 110,0.,1100.);
  h_tickmark_height_3 = fs->make<TH1F>("tickmark_height_3", "Tickmark height with AOH gain at 3", 110,0.,1100.);

  h_tickmark_replaced_0 = fs->make<TH1F>("tickmark_replaced_0", "Tickmark height with AOH gain at 0", 110,0.,1100.);
  h_tickmark_replaced_1 = fs->make<TH1F>("tickmark_replaced_1", "Tickmark height with AOH gain at 1", 110,0.,1100.);
  h_tickmark_replaced_2 = fs->make<TH1F>("tickmark_replaced_2", "Tickmark height with AOH gain at 2", 110,0.,1100.);
  h_tickmark_replaced_3 = fs->make<TH1F>("tickmark_replaced_3", "Tickmark height with AOH gain at 3", 110,0.,1100.);

  h_tickmark_recovered_0 = fs->make<TH1F>("tickmark_recovered_0", "Tickmark height with AOH gain at 0", 110,0.,1100.);
  h_tickmark_recovered_1 = fs->make<TH1F>("tickmark_recovered_1", "Tickmark height with AOH gain at 1", 110,0.,1100.);
  h_tickmark_recovered_2 = fs->make<TH1F>("tickmark_recovered_2", "Tickmark height with AOH gain at 2", 110,0.,1100.);
  h_tickmark_recovered_3 = fs->make<TH1F>("tickmark_recovered_3", "Tickmark height with AOH gain at 3", 110,0.,1100.);

  h_height_ratio_at_DG0  = fs->make<TH1F>("height_ratio_at_DG0",  "Tickmark height ratio with #Delta AOHGain==0",400,0.,+2.);
  h_height_ratio_at_DGp1 = fs->make<TH1F>("height_ratio_at_DGp1", "Tickmark height ratio with #Delta AOHGain==1",400,0.,+2.);
  h_height_ratio_at_DGp2 = fs->make<TH1F>("height_ratio_at_DGp2", "Tickmark height ratio with #Delta AOHGain==2",400,0.,+2.);
  h_height_ratio_at_DGp3 = fs->make<TH1F>("height_ratio_at_DGp3", "Tickmark height ratio with #Delta AOHGain==3",400,0.,+2.);
  h_height_ratio_at_DGm1 = fs->make<TH1F>("height_ratio_at_DGm1", "Tickmark height ratio with #Delta AOHGain==1",400,0.,+2.); 
  h_height_ratio_at_DGm2 = fs->make<TH1F>("height_ratio_at_DGm2", "Tickmark height ratio with #Delta AOHGain==2",400,0.,+2.);
  h_height_ratio_at_DGm3 = fs->make<TH1F>("height_ratio_at_DGm3", "Tickmark height ratio with #Delta AOHGain==3",400,0.,+2.);

  for (unsigned int i=0; i<15; i++) {
    Summary s;
    s.det_type = i;
    std::string tag        = "noise_vs_gain_"+s.detType();
    std::string tit        = "Noise versus Gain for "+s.detType();
    h_noise_vs_gain[i]     = fs->make<TH2F>(tag.c_str(),tit.c_str(),100,0.,2.,124,0.,31);
    std::string btag       = "bad_noise_vs_gain_"+s.detType();
    std::string btit       = "Noise versus Gain of bad channels for "+s.detType();
    h_bad_channels[i]      = fs->make<TH2F>(btag.c_str(),btit.c_str(),100,0.,2.,124,0.,31);
    std::string stag       = "surveyed_noise_vs_gain_"+s.detType();
    std::string stit       = "Noise versus Gain of surveyed channels for "+s.detType();
    h_surveyed_channels[i] = fs->make<TH2F>(stag.c_str(),stit.c_str(),100,0.,2.,124,0.,31);
    std::string ntag       = "noisy_noise_vs_gain_"+s.detType();
    std::string ntit       = "Noise versus Gain of noisy channels for "+s.detType();
    h_noisy_channels[i]    = fs->make<TH2F>(ntag.c_str(),ntit.c_str(),100,0.,2.,124,0.,31);
    std::string gtag       = "good_noise_vs_gain_"+s.detType();
    std::string gtit       = "Noise versus Gain of channels from "+s.detType();
    h_good_channels[i]     = fs->make<TH2F>(gtag.c_str(),gtit.c_str(),100,0.,2.,124,0.,31);
    std::string etag       = "extreme_noise_vs_gain_"+s.detType();
    std::string etit       = "Noise versus Gain of extreme channels from "+s.detType();
    h_badAOHgs_channels[i] = fs->make<TH2F>(etag.c_str(),etit.c_str(),100,0.,2.,124,0.,31);
  }
}

void
SiStripApvGainFromFileBuilder::endJob() {
  // Register to the TFileService
  edm::Service<TFileService> fs;

  //fitting noise histograms
  for(unsigned int i=1; i<15;i++) {
    TProfile* p = fs->make<TProfile>( *h_noise_vs_gain[i]->ProfileX() );
    p->Fit( "pol1","","",0.85,1.3 );
    TF1* f = p->GetFunction("pol1");
    float a     = f->GetParameter(1);
    //float err_a = f->GetParError(1)*5.;
    float b     = f->GetParameter(0);
    //float err_b = f->GetParError(0)*5.;

    //checking the status of all the "bad" channels
    std::vector<Summary>::iterator s = summary_.begin();
    while ( s!=summary_.end() ) {
      Summary& summary = (*s);

      float noiseFromTiming    = a*summary.gain_in_db +b;
      float noiseFromReference = a*(summary.reference_height/640.) +b;

      float noiseRatio = summary.noise/noiseFromTiming;

      bool correctDetType = ( summary.det_type==i );

      //bool extremeHeight = ( summary.tickmark_height>830. || summary.tickmark_height<545. );

      if ( correctDetType && summary.tickmark_status=='L' && (noiseFromTiming==0. ||
                                            std::fabs((noiseFromTiming-summary.noise)/noiseFromTiming)>0.1) ) {
        // Set the status as BAD and move the summary to ex_summary list
        summary.tickmark_status='B';
        Summary new_summary = summary;
        s = summary_.erase( s );

        ex_summary_.push_back( new_summary );

        h_bad_channels[i]->Fill( new_summary.gain_in_db, new_summary.noise );

      } else if ( correctDetType && (std::fabs(noiseFromTiming)-summary.noise)/std::fabs(noiseFromTiming)>0.35) {
          // set the status as Noisy
          summary.tickmark_status='N';
          h_noisy_channels[i]->Fill( summary.gain_in_db, summary.noise );
          if ( (std::fabs(noiseFromReference)-summary.noise)/std::fabs(noiseFromReference)>0.35 ) {
            summary.status_wrt_previous_run = '=';
          } else {
            summary.status_wrt_previous_run = '-';
          }
          s++;
      } else if ( correctDetType && summary.tickmark_status=='V' && std::fabs( 1.-noiseRatio )>0.05 ) {
        // nothing to do, discrepancy with reference is more than 5%, keep eyes open
        h_surveyed_channels[i]->Fill( summary.gain_in_db, summary.noise );
        s++;
      } else if ( correctDetType && (!is_aoh_gain_ok(summary.tickmark_height,summary.aoh_gain)) ) {
        // this one has wrong AOH gain setting
        summary.tickmark_status='E';
        h_badAOHgs_channels[i]->Fill( summary.gain_in_db, summary.noise );
        if ( !is_aoh_gain_ok(summary.reference_height,summary.reference_aoh_gain) ) summary.status_wrt_previous_run='=';
        else summary.status_wrt_previous_run='-';
        s++;
      }  else {
        if ( correctDetType) {
          summary.tickmark_status='G';
          h_good_channels[i]->Fill( summary.gain_in_db, summary.noise );
        }
        s++;
      }

    }
  }

  if( outputSummary_ ) output_summary();

  try {
    this->output_ratio_maps();
    this->output_topology_maps();
  } catch ( std::exception& e ) {
    std::cerr << e.what() << std::endl;
  }


}


void SiStripApvGainFromFileBuilder::analyze(const edm::Event& evt, const edm::EventSetup& iSetup){

  //unsigned int run=evt.id().run();

  edm::LogInfo("Workflow") << "@SUB=analyze" <<  "Insert SiStripApvGain Data." << std::endl;
  this->read_tickmark();
  this->read_summary(recovery_);
  if ( !this->read_noise_from_root( nfp_.fullPath().c_str(), &noise_ ) ) exit(1);

  if (outputMaps_) {
    try {
      this->output_tickmark_maps(&gains_,"tickmark_heights");
      this->output_tickmark_maps(&negative_gains_,"negative_tickmark");
      this->output_tickmark_maps(&null_gains_,"zero_tickmark");
    } catch ( std::exception& e ) {
      std::cerr << e.what() << std::endl;
    }
  }

  // Retrieve the SiStripDetCabling description
  iSetup.get<SiStripDetCablingRcd>().get( detCabling_ );

  //Retrieve the RunInfo object
  edm::ESHandle<RunInfo> runInfo;
  iSetup.get<RunInfoRcd>().get( runInfo );

  siStripQuality_->add( detCabling_.product() );
  siStripQuality_->add( runInfo.product() );
  siStripQuality_->cleanUp();
  siStripQuality_->fillBadComponents();

  // APV gain record to be filled with gains and delivered into the database.
  SiStripApvGain* obj = new SiStripApvGain();

  SiStripDetInfoFileReader reader(gfp_.fullPath());
  
  const std::map<uint32_t, SiStripDetInfoFileReader::DetInfo > DetInfos  = reader.getAllData();

  LogTrace("Debug") << "  det id  |APVOF| CON |APVON| FED |FEDCH|i2cAd|tickHeig|AOH gain|AOH bias|" << std::endl;

  for(std::map<uint32_t, SiStripDetInfoFileReader::DetInfo >::const_iterator it = DetInfos.begin(); it != DetInfos.end(); it++){    
    // check if det id is correct and if it is actually cabled in the detector
    if( it->first==0 || it->first==0xFFFFFFFF ) {
      edm::LogError("DetIdNotGood") << "@SUB=analyze" << "Wrong det id: " << it->first 
                                    << "  ... neglecting!" << std::endl;
      continue;
    }


    // For the cabled det_id retrieve the number of APV connected
    // to the module with the FED cabling
    uint16_t nAPVs = 0;
    const std::vector<const FedChannelConnection*> connection = detCabling_->getConnections(it->first);
    for (unsigned int ca = 0; ca<connection.size(); ca++) {
      if ( connection[ca]!=0 )  {
        nAPVs+=( connection[ca] )->nApvs();
        break;
      }
    }

    // check consistency among FED cabling and ideal cabling, exit on error
    if (connection.size()!=0 && nAPVs != (uint16_t)it->second.nApvs) {
      edm::LogError("SiStripCablingError") << "@SUB=analyze"
                     << "det id " << it->first << ": APV number from FedCabling (" << nAPVs
                     << ") is different from the APV number retrieved from the ideal cabling ("
                     << it->second.nApvs << ")." << std::endl;
      throw("Inconsistency on the number of APVs.");
    }


    // eventually separate the processing for the module that are fully
    // uncabled. This is worth only if we decide not tu put the record
    // in the DB for the uncabled det_id.
    //if( !detCabling_->IsConnected(it->first) ) {
    //
    //  continue;
    //}
 

    //Gather the APV online id
    //std::vector<std::pair<int,float>> tickmark_for_detId(it->second.nApvs,std::pair<int,float>(-1,999999.));
    std::vector<std::pair<int,tickmark_info>> 
        tickmark_for_detId(it->second.nApvs,std::pair<int,tickmark_info>(-1,tickmark_info() ));
    std::vector<std::pair<int,tickmark_info>>
        tickmark_reference(it->second.nApvs,std::pair<int,tickmark_info>(-1,tickmark_info() ));
    std::vector<std::pair<int,noise_info>>
        noise_for_detId(it->second.nApvs,std::pair<int,noise_info>(-1,noise_info() ));

    for (unsigned int ca = 0; ca<connection.size(); ca++) {
      if( connection[ca]!=0 ) {
        uint16_t id1 = (connection[ca])->i2cAddr( 0 )%32;
        uint16_t id2 = (connection[ca])->i2cAddr( 1 )%32;
        tickmark_for_detId[ online2offline(id1,it->second.nApvs) ].first = id1;
        tickmark_for_detId[ online2offline(id2,it->second.nApvs) ].first = id2;
        noise_for_detId[ online2offline(id1,it->second.nApvs) ].first = id1;
        noise_for_detId[ online2offline(id2,it->second.nApvs) ].first = id2;
      }
    }

    height_from_maps(it->first, it->second.nApvs, tickmark_for_detId);
    height_from_reference(it->first, it->second.nApvs, tickmark_reference);
    noise_from_maps(it->first, it->second.nApvs, noise_for_detId);

    std::vector<float> theSiStripVector;

    // Fill the gain in the DB object, apply the logic for the dummy values
    for(unsigned short j=0; j<it->second.nApvs; j++){
      Summary summary;
      summary.det_id = it->first;
      summary.det_type = SiStripDetId(it->first).moduleGeometry();
      summary.offlineAPV_id = j;
      summary.onlineAPV_id  = tickmark_for_detId.at(j).first;
      summary.is_connected = false;
      summary.strip_length = reader.getNumberOfApvsAndStripLength(it->first).second;
      summary.noise    = noise_for_detId.at(j).second.noise;
      summary.pedestal = noise_for_detId.at(j).second.pedestal;  
      summary.FED_id = -1;
      summary.FED_ch = -1;
      summary.i2cAdd = -1;
      summary.FEC_crate = -1;
      summary.FEC_slot  = -1;
      summary.FEC_ring  = -1;
      summary.CCU_addr  = -1;
      summary.CCU_chan  = -1;

      summary.channel_status = ( siStripQuality_->IsApvBad(summary.det_id,summary.offlineAPV_id) )? 'd' : 'e';

      for (unsigned int ca = 0; ca<connection.size(); ca++) { 
        if( connection[ca]!=0 && (connection[ca])->i2cAddr( j%2 )%32==summary.onlineAPV_id ) {
          summary.is_connected = (connection[ca])->isConnected();
          summary.FED_id = (connection[ca])->fedId();
          summary.FED_ch = (connection[ca])->fedCh();
          summary.i2cAdd = (connection[ca])->i2cAddr( j%2 );
          summary.FEC_crate = (connection[ca])->fecCrate();
          summary.FEC_slot  = (connection[ca])->fecSlot();
          summary.FEC_ring  = (connection[ca])->fecRing();
          summary.CCU_addr  = (connection[ca])->ccuAddr();
          summary.CCU_chan  = (connection[ca])->ccuChan();
        }
      }

      try {
        float height           = tickmark_for_detId[j].second.tickmark_height;
        float reference_height = tickmark_reference[j].second.tickmark_height;

        summary.tickmark_height = height;
        summary.reference_height = reference_height;
        summary.aoh_gain = tickmark_for_detId[j].second.aoh_gain;
        summary.aoh_bias = tickmark_for_detId[j].second.aoh_bias;
        summary.reference_aoh_gain = tickmark_reference[j].second.aoh_gain;
        summary.reference_aoh_bias = tickmark_reference[j].second.aoh_bias;
        LogTrace("Debug") << it->first << "  " << std::setw(3) << j << "   "
                          << std::setw(3) << connection.size() << "   "
                          << std::setw(3) << summary.onlineAPV_id << "    "
                          << std::setw(3) << summary.FED_id << "   "
                          << std::setw(3) << summary.FED_ch << "   "
                          << std::setw(3) << summary.i2cAdd << "   "
                          << std::setw(7) << summary.tickmark_height 
                          << std::setw(7) << summary.aoh_gain
                          << std::setw(7) << summary.aoh_bias << std::endl;

        summary.gain_in_db = 1.;        // use 1. as default value to avoid possible crashes at reconstruction level
        summary.tickmark_status = 'G';  // use good as default
        summary.recovery_status = 'U';  // unrecovered gain
        if ( is_usable(summary.reference_height) )  summary.status_wrt_previous_run = '=';
        else                                        summary.status_wrt_previous_run = '+';

        set_tickmark_status(summary);

        if ( height!=999999. ) {
          summary.is_scanned = true;

          if ( summary.aoh_gain==0 ) h_tickmark_height_0->Fill(height,1.);
          else if ( summary.aoh_gain==1 ) h_tickmark_height_1->Fill(height,1.);
          else if ( summary.aoh_gain==2 ) h_tickmark_height_2->Fill(height,1.);
          else if ( summary.aoh_gain==3 ) h_tickmark_height_3->Fill(height,1.);


          if( height>(float)heightThreshold_ ) {
            set_gain(summary, height, 'U');

            float recovery_height = apply_recovery(summary);
            if( doRecovery_ && recovery_height!=0. ) set_gain(summary, recovery_height, 'R');

            if ( !summary.is_connected ) ex_summary_.push_back( summary );
            else summary_.push_back( summary );

            // Fill histogram to monitor the quality
            h_noise_vs_gain[summary.det_type]->Fill(height/640.,summary.noise);

            if ( reference_height!=999999. && reference_height>0. ) {
              float ratio = height / reference_height;
              int dAOHGain = summary.aoh_gain - tickmark_reference[j].second.aoh_gain;
              if      (dAOHGain==0)  h_height_ratio_at_DG0->Fill(ratio,1.0);
              else if (dAOHGain==1)  h_height_ratio_at_DGp1->Fill(ratio,1.0);
              else if (dAOHGain==2)  h_height_ratio_at_DGp2->Fill(ratio,1.0);
              else if (dAOHGain==3)  h_height_ratio_at_DGp3->Fill(ratio,1.0);
              else if (dAOHGain==-1) h_height_ratio_at_DGm1->Fill(ratio,1.0);
              else if (dAOHGain==-2) h_height_ratio_at_DGm2->Fill(ratio,1.0);
              else if (dAOHGain==-3) h_height_ratio_at_DGm3->Fill(ratio,1.0);

            }

          } else {
            float recovery_height = apply_recovery(summary);
            if ( height==0. ) {
              // channels that were off during the timing scan are processed here
              if( doRecovery_ && recovery_height!=0. ) set_gain(summary,recovery_height,'R');
              else if( putDummyIntoOffChannels_ )      set_gain(summary,dummyAPVGain_,'D');
              ex_summary_.push_back( summary );
            }
            else if (height<0. )   {
              // channels sending corrupted data during the timing scan are processed here
              if( doRecovery_ && recovery_height!=0. ) set_gain(summary,recovery_height,'R');
              else if( putDummyIntoBadChannels_ )      set_gain(summary,dummyAPVGain_,'D');
              ex_summary_.push_back( summary );
            }
            else {
              h_noise_vs_gain[summary.det_type]->Fill(height/640.,summary.noise);
              set_gain(summary,height,'U');
              summary.tickmark_status = 'L';
              summary.status_wrt_previous_run = ( is_usable( summary.reference_height ) )? '-' : '=';

              // channles whose tickmark height is below threshold are processed here
              if( doRecovery_ && recovery_height!=0. ) set_gain(summary,recovery_height,'R');
              else if ( putDummyIntoLowChannels_ )     set_gain(summary,dummyAPVGain_,'D');
              summary_.push_back( summary );
            }
          } 
        } else {
          // channels which are not scanned during the timing run are processed here
          summary.is_scanned = false;
          float recovery_height = apply_recovery(summary);
          if( !summary.is_connected ) {
            if( doRecovery_ && recovery_height!=0. ) set_gain(summary,recovery_height,'R');
            else if( putDummyIntoUncabled_ )         set_gain(summary,dummyAPVGain_,'D');
          } else {
            if( doRecovery_ && recovery_height!=0. ) set_gain(summary,recovery_height,'R');
            else if( putDummyIntoUnscanned_ )        set_gain(summary,dummyAPVGain_,'D');
          }
          ex_summary_.push_back( summary );
        }

        theSiStripVector.push_back( summary.gain_in_db );

      } catch ( std::exception& e ) {
        std::cerr << e.what() << std::endl;
        edm::LogError("MappingError") << "@SUB=analyze" << "Job end prematurely." << std::endl;
        exit(1);
      }
    }
  	    
    SiStripApvGain::Range range(theSiStripVector.begin(),theSiStripVector.end());
    if ( ! obj->put(it->first,range) )
      edm::LogError("IndexError")<< "@SUB=analyze" << "detid already exists." << std::endl;
  }


  
  //End now write sistrip gain data in DB
  edm::Service<cond::service::PoolDBOutputService> mydbservice;
  
  if( mydbservice.isAvailable() ){
    if( mydbservice->isNewTagRequest("SiStripApvGainRcd") ){
      mydbservice->createNewIOV<SiStripApvGain>(obj,mydbservice->beginOfTime(),mydbservice->endOfTime(),"SiStripApvGainRcd");      
    } else {
      mydbservice->appendSinceTime<SiStripApvGain>(obj,mydbservice->currentTime(),"SiStripApvGainRcd");      
    }
  }else{
    edm::LogError("DBServiceNotAvailable") << "@SUB=analyze" << "DB Service is unavailable" << std::endl;
  }
}

bool
SiStripApvGainFromFileBuilder::read_summary(std::vector<Summary>& list) {
  list.clear();
  const char* filename = recoveryList_.fullPath().c_str();

  std::ifstream summary_list ( filename );
  if ( !summary_list.is_open() ) {
    edm::LogError("FileNotFound") << "@SUB=read_summary" << "File with summary list "
                                  << filename << " cannot be opened!" << std::endl;
    return false;
  }

    std::string header = "  det id |APV| Connected |FED|CH|FEC|FECs|FECr|CCU|CCUc|i2cA|APVON|  noise  |tickHeig|AOH|Bias|DBgain|refHeigh|AOH|Bias|Status|";
    std::string ruler  = "---------+---+-----------+---+--+---+----+----+---+----+----+-----+---------+--------+---+----+------+--------+---+----+------+";

  int count = -1;

  for ( ; ; ) {
    count++;
    Summary s;
    std::string con1,con2;
    //std::string scanning;

    summary_list >> s.det_id          >> s.offlineAPV_id   >> con1 >> con2       >> s.FED_id             >> s.FED_ch 
                 >> s.FEC_crate       >> s.FEC_slot        >> s.FEC_ring         >> s.CCU_addr           >> s.CCU_chan
                 >> s.i2cAdd          >> s.onlineAPV_id    >> s.noise            >> s.tickmark_height    >> s.aoh_gain
                 >> s.aoh_bias        >> s.gain_in_db      >> s.reference_height >> s.reference_aoh_gain >> s.reference_aoh_bias
                 >> s.tickmark_status >> s.channel_status  >> s.recovery_status  >> s.status_wrt_previous_run;

    s.is_connected = (con1.compare("NOT")==0)? false : true;
    s.det_type     = s.detType( con2.c_str() );
    s.is_scanned   = (s.tickmark_height!=999999.)?   true : false;

    if ( ! (summary_list.eof() || summary_list.fail()) ) {

      if (count==0) {
        edm::LogVerbatim("Workflow") << "Reading " << filename << " for gathering the summary list" << std::endl;
        edm::LogVerbatim("Workflow") << header << std::endl;
        edm::LogVerbatim("Workflow") << ruler << std::endl;
      }
      std::stringstream line;
      format_summary(line, s);

      edm::LogVerbatim("Workflow") << line.str() << std::endl;

      // insert the gain value in the map
      recovery_.push_back( s );

    } else if (summary_list.eof()) {
      edm::LogInfo("Workflow") << "@SUB=read_summary" << "EOF of " << filename << " reached." << std::endl;
      break;
    } else if (summary_list.fail()) {
      edm::LogError("FileiReadError") << "@SUB=read_summary" << "error while reading " << filename << std::endl;
      break;
    }
  }

  summary_list.close();
  return true;
}
     
void
SiStripApvGainFromFileBuilder::read_tickmark() {

  // clear internal maps
  for_each(gains_.begin(), gains_.end(), CleanUp);
  for_each(negative_gains_.begin(), negative_gains_.end(), CleanUp);
  for_each(null_gains_.begin(), null_gains_.end(), CleanUp);

  for_each(ref_gains_.begin(), ref_gains_.end(), CleanUp);


  // Connect file for input
  const char* filename = tfp_.fullPath().c_str();
  bool read_ok = false;

  // check if a root file is connected
  if ( ends_with(filename,".root") ) {
    read_ok = this->read_tickmark_from_root(filename,&gains_,&negative_gains_,&null_gains_);
  } else {
    read_ok = this->read_tickmark_from_ascii(filename,&gains_,&negative_gains_,&null_gains_);
  }

  if (read_ok==false) {
    exit(1);
  }


  // Connect the reference file if any
  filename = rfp_.fullPath().c_str(); 
  if ( ends_with(filename,".root") ) {
    read_ok = this->read_tickmark_from_root(filename,&ref_gains_,&ref_gains_,&ref_gains_);
  } else {
    read_ok = this->read_tickmark_from_ascii(filename,&ref_gains_,&ref_gains_,&ref_gains_);
  }

  if (read_ok==false) {
    exit(1);
  }
}


bool
SiStripApvGainFromFileBuilder::read_tickmark_from_root(const char* filename,
                                                       std::vector<Gain*>* gains,
                                                       std::vector<Gain*>* negative_gains,
                                                       std::vector<Gain*>* null_gains) {

    TFile* tickmark_heights = TFile::Open( filename );
    if ( !tickmark_heights->IsOpen() ) {
      edm::LogError("FileNotFound") << "@SUB=read_tickmark_from_root" << "File with thickmark height "
                                    << filename << " cannot be opened!" << std::endl;
      return false;
    }

    // gather the ntuple
    TTree* ntuple = (TTree*) tickmark_heights->Get("tree");
    if (ntuple==0) {
      edm::LogError("TreeNotFound") << "@SUB=read_tickmark_from_root" 
                                    << "TTree with tickmark scan data not found in " << filename << std::endl;
      return false;
    }

    //connect branches
    double det_id     = 999999.;
    double apv_id     = 999999.;
    double apv_height = 999999.;
    double aoh_gain   = 256.;
    double aoh_bias   = 256.;
    if ( ntuple->GetBranch("DETID")!=0 ) ntuple->SetBranchAddress("DETID",&det_id);
    if ( ntuple->GetBranch("APVID")!=0 ) ntuple->SetBranchAddress("APVID",&apv_id);
    if ( ntuple->GetBranch("HEIGHT")!=0 ) ntuple->SetBranchAddress("HEIGHT",&apv_height);
    if ( ntuple->GetBranch("GAIN")!=0 ) ntuple->SetBranchAddress("GAIN",&aoh_gain);
    if ( ntuple->GetBranch("BIAS")!=0 ) ntuple->SetBranchAddress("BIAS",&aoh_bias);

    int entries = ntuple->GetEntries();

    int count = -1;

    for( int i=0; i<entries; i++) {
      ntuple->GetEntry(i);
      count++;

      if (count==0) {
        LogTrace("Debug") << "Reading " << filename << " for gathering the tickmark heights" << std::endl;
        LogTrace("Debug") << "|  Det Id   |  APV Id  |   Tickmark   |   AOH GAIN   |   AOH BIAS   " << std::endl;
        LogTrace("Debug") << "+-----------+----------+--------------+--------------+--------------" << std::endl;
      }
      LogTrace("Debug") << std::setw(11) << det_id
                        << std::setw(8)  << apv_id
                        << std::setw(14) << apv_height
                        << std::setw(14) << aoh_gain
                        << std::setw(14) << aoh_bias << std::endl;

      // retrieve the map corresponding to the gain collection
      Gain* map = 0;
      if ( apv_height>0. ) {
        map = get_gain_map(gains, apv_id);
      } else if ( apv_height<0. ) {
        map = get_gain_map(negative_gains, apv_id);
      } else if ( apv_height==0.) {
        map = get_gain_map(null_gains, apv_id);
      }

      // insert the gain value in the map
      tickmark_info tickmarkInfo(apv_height, aoh_gain, aoh_bias);
      
      std::pair<Gain::iterator,bool> ret = map->insert( std::pair<uint32_t,tickmark_info>(det_id, tickmarkInfo) );
      if ( ret.second == false ) {
        edm::LogError("MapError") << "@SUB=read_tickmark_from_root" << "Cannot not insert gain for detector id "
                                  << det_id << " into the internal map: detector id already in the map." << std::endl;
      }

    }

    return true;
}


bool
SiStripApvGainFromFileBuilder::read_tickmark_from_ascii(const char* filename,
                                                       std::vector<Gain*>* gains,
                                                       std::vector<Gain*>* negative_gains,
                                                       std::vector<Gain*>* null_gains) {

  std::ifstream thickmark_heights ( filename );
  if ( !thickmark_heights.is_open() ) {
    edm::LogError("FileNotFound") << "@SUB=read_tickmark_from_ascii" << "File with thickmark height " 
                                  << filename << " cannot be opened!" << std::endl;
    return false;
  }


  // read file and fill internal map
  uint32_t det_id = 0;
  uint32_t APV_id = 0;
  float    tick_h = 0.;

  int count = -1;

  for ( ; ; ) {
    count++;
    thickmark_heights >> det_id >> APV_id >> tick_h;

    if ( ! (thickmark_heights.eof() || thickmark_heights.fail()) ) {

      if (count==0) {
        LogTrace("Debug") << "Reading " << filename << " for gathering the tickmark heights" << std::endl;
        LogTrace("Debug") << "|  Det Id   |  APV Id  |  Tickmark" << std::endl;
        LogTrace("Debug") << "+-----------+----------+----------" << std::endl;
      }
      LogTrace("Debug") << std::setw(11) << det_id
                        << std::setw(8)  << APV_id 
                        << std::setw(14) << tick_h << std::endl;

      // retrieve the map corresponding to the gain collection
      Gain* map = 0;
      if ( tick_h>0. ) {
        map = get_gain_map(gains, APV_id); 
      } else if ( tick_h<0. ) {
        map = get_gain_map(negative_gains, APV_id);
      } else if ( tick_h==0.) {
        map = get_gain_map(null_gains, APV_id);
      }

      // insert the gain value in the map
      std::pair<Gain::iterator,bool> ret = map->insert( std::pair<uint32_t,tickmark_info>(det_id,tickmark_info(tick_h) ));
      if ( ret.second == false ) {
        edm::LogError("MapError") << "@SUB=read_tickmark_from_ascii" << "Cannot not insert gain for detector id "
                                  << det_id << " into the internal map: detector id already in the map." << std::endl;
      }

    } else if (thickmark_heights.eof()) {
      edm::LogInfo("Workflow") << "@SUB=read_tickmark_from_ascii" << "EOF of " << filename << " reached." << std::endl;
      break;
    } else if (thickmark_heights.fail()) {
      edm::LogError("FileiReadError") << "@SUB=read_tickmark_from_ascii" << "error while reading " << filename << std::endl;
      break;
    }
  }

  thickmark_heights.close();
  return true;
}

bool
SiStripApvGainFromFileBuilder::read_noise_from_root( const char* filename, std::vector<Noise*>* noise) {

  TFile* noise_and_pedestal = TFile::Open( filename );
  if ( !noise_and_pedestal->IsOpen() ) {
    edm::LogError("FileNotFound") << "@SUB=read_noise_from_root" << "File with noise "
                                  << filename << " cannot be opened!" << std::endl;
    return false;
  }

  // gather the ntuple
  TTree* ntuple = (TTree*) noise_and_pedestal->Get("tree");
  if (ntuple==0) {
    edm::LogError("TreeNotFound") << "@SUB=read_noise_from_root"
                                  << "TTree with tickmark scan data not found in " << filename << std::endl;
    return false;
  }

  //connect branches
  double det_id       = 999999.;
  double apv_id       = 999999.;
  double apv_noise    = 999999.;
  double apv_pedestal = 256.;
  if ( ntuple->GetBranch("DETID")!=0 ) ntuple->SetBranchAddress("DETID",&det_id);
  if ( ntuple->GetBranch("APVID")!=0 ) ntuple->SetBranchAddress("APVID",&apv_id);
  if ( ntuple->GetBranch("NOISEMEAN")!=0 ) ntuple->SetBranchAddress("NOISEMEAN",&apv_noise);
  if ( ntuple->GetBranch("PEDSMEAN")!=0 ) ntuple->SetBranchAddress("PEDSMEAN",&apv_pedestal);

  int entries = ntuple->GetEntries();

  int count = -1;

  for( int i=0; i<entries; i++) {
    ntuple->GetEntry(i);
    count++;

    if (count==0) {
      LogTrace("Debug") << "Reading " << filename << " for gathering the tickmark heights" << std::endl;
      LogTrace("Debug") << "|  Det Id   |  APV Id  |     Noise    |   Pedestal   " << std::endl;
      LogTrace("Debug") << "+-----------+----------+--------------+--------------" << std::endl;
    }
    LogTrace("Debug") << std::setw(11) << det_id
                      << std::setw(8)  << apv_id
                      << std::setw(14) << apv_noise
                      << std::setw(14) << apv_pedestal << std::endl;

    // retrieve the map corresponding to the gain collection
    Noise* map = 0;
    map = get_noise_map(noise, apv_id);

    // insert the gain value in the map
    noise_info noiseInfo(apv_noise, apv_pedestal);

    std::pair<Noise::iterator,bool> ret = map->insert( std::pair<uint32_t,noise_info>(det_id, noiseInfo) );
    if ( ret.second == false ) {
      edm::LogError("MapError") << "@SUB=read_tickmark_from_root" << "Cannot not insert gain for detector id "
                                << det_id << " into the internal map: detector id already in the map." << std::endl;
    }

  }

  return true;
}


SiStripApvGainFromFileBuilder::Gain* 
SiStripApvGainFromFileBuilder::get_gain_map(std::vector<Gain*>* maps, int onlineAPV_id) {
    Gain* map=0;
    if( onlineAPV_id<0 || onlineAPV_id>5 ) return map;

    try {
        map = maps->at(onlineAPV_id);
    } catch (std::out_of_range) {
        if ( maps->size()<static_cast<unsigned int>(onlineAPV_id) ) maps->resize( onlineAPV_id );
        maps->insert( maps->begin()+onlineAPV_id, new Gain() );
        map = (*maps)[onlineAPV_id];
    }

    if ( map==0 ) {
        (*maps)[onlineAPV_id] = new Gain();
        map = (*maps)[onlineAPV_id];
    }

    return map;
}


SiStripApvGainFromFileBuilder::Noise*
SiStripApvGainFromFileBuilder::get_noise_map(std::vector<Noise*>* maps, int onlineAPV_id) {
  Noise* map=0;
  if( onlineAPV_id<0 || onlineAPV_id>5 ) return map;

  try {
      map = maps->at(onlineAPV_id);
  } catch (std::out_of_range) {
      if ( maps->size()<static_cast<unsigned int>(onlineAPV_id) ) maps->resize( onlineAPV_id );
      maps->insert( maps->begin()+onlineAPV_id, new Noise() );
      map = (*maps)[onlineAPV_id];
  }

  if ( map==0 ) {
      (*maps)[onlineAPV_id] = new Noise();
      map = (*maps)[onlineAPV_id];
  }

  return map;
}


void SiStripApvGainFromFileBuilder::output_topology_maps() const {
    // check if directory exist, if not create it
    struct stat sb;
    if ( !(stat("tracker_maps",&sb)==0 && S_ISDIR(sb.st_mode)) ) {
      if( mkdir("tracker_maps",0777)==-1 ) throw "cannot create the tracker_maps directory!";
    }

    // open output file
    std::ofstream* oNSY[3] = {0, 0, 0};
    std::set< std::pair<uint32_t, float> > NSY[3];

    std::ofstream* oEXT[3] = {0, 0, 0};
    std::set< std::pair<uint32_t, float> > EXT[3];

    std::ofstream* oTBC[3] = {0, 0, 0};
    std::set< std::pair<uint32_t, float> > TBC[3];

    std::ofstream* oBAD[3] = {0, 0, 0};
    std::set< std::pair<uint32_t, float> > BAD[3];

    std::ofstream* oUNC[3] = {0, 0, 0};
    std::set< std::pair<uint32_t, float> > UNC[3];


    for (unsigned int i=0; i<3; i++) {
      std::stringstream name_NSY,err_NSY;
      name_NSY << "tracker_maps/Noisy_APVpair" << i << ".txt";
      err_NSY  << "cannot open file " << name_NSY.str();
      oNSY[i] = new std::ofstream(name_NSY.str(), std::ofstream::trunc);
      if ( !oNSY[i]->is_open() ) throw  err_NSY.str();

      std::stringstream name_EXT,err_EXT;
      name_EXT << "tracker_maps/WrongAOH_APVpair" << i << ".txt";
      err_EXT  << "cannot open file " << name_EXT.str();
      oEXT[i] = new std::ofstream(name_EXT.str(), std::ofstream::trunc);
      if ( !oEXT[i]->is_open() ) throw  err_EXT.str();

      std::stringstream name_TBC,err_TBC;
      name_TBC << "tracker_maps/ToBeChecked_APVpair" << i << ".txt";
      err_TBC  << "cannot open file " << name_TBC.str();
      oTBC[i] = new std::ofstream(name_TBC.str(), std::ofstream::trunc);
      if ( !oTBC[i]->is_open() ) throw  err_TBC.str();

      std::stringstream name_BAD,err_BAD;
      name_BAD << "tracker_maps/Bad_APVpair" << i << ".txt";
      err_BAD  << "cannot open file " << name_BAD.str();
      oBAD[i] = new std::ofstream(name_BAD.str(), std::ofstream::trunc);
      if ( !oBAD[i]->is_open() ) throw  err_BAD.str();

      std::stringstream name_UNC,err_UNC;
      name_UNC << "tracker_maps/UncabledAndScanned_APVpair" << i << ".txt";
      err_UNC  << "cannot open file " << name_UNC.str();
      oUNC[i] = new std::ofstream(name_UNC.str(), std::ofstream::trunc);
      if ( !oUNC[i]->is_open() ) throw  err_UNC.str();
    }

    for(unsigned int i=0; i<summary_.size(); i++) {
      Summary s = summary_[i];
      int APV_pairIndex = s.offlineAPV_id/2;
      std::pair<uint32_t,float> p = std::pair<uint32_t,float>(s.det_id,s.tickmark_height);
      if ( (!s.is_connected) && (s.is_scanned) ) UNC[APV_pairIndex].insert( p );
      else if ( s.tickmark_status=='B' )         BAD[APV_pairIndex].insert( p );
      else if ( s.tickmark_status=='V' )         TBC[APV_pairIndex].insert( p );
      else if ( s.tickmark_status=='E' )         EXT[APV_pairIndex].insert( p );
      else if ( s.tickmark_status=='N' )         NSY[APV_pairIndex].insert( p );
    }

    for(unsigned int i=0; i<ex_summary_.size(); i++) {
      Summary s = ex_summary_[i];
      int APV_pairIndex = s.offlineAPV_id/2;
      std::pair<uint32_t,float> p = std::pair<uint32_t,float>(s.det_id,s.tickmark_height);
      if ( (!s.is_connected) && (s.is_scanned) ) UNC[APV_pairIndex].insert( p );
      else if ( s.tickmark_status=='B' )         BAD[APV_pairIndex].insert( p );
      else if ( s.tickmark_status=='V' )         TBC[APV_pairIndex].insert( p );
      else if ( s.tickmark_status=='E' )         EXT[APV_pairIndex].insert( p );
      else if ( s.tickmark_status=='N' )         NSY[APV_pairIndex].insert( p );
    }

    // dump the set content into the output files
    for(unsigned int APV=0; APV<3; APV++) {
      std::set< std::pair<uint32_t,float>> uncabled = UNC[APV];
      std::set< std::pair<uint32_t,float>>::iterator it;
      for(it=uncabled.begin();it!=uncabled.end();it++) (*oUNC[APV]) << (*it).first << "   " << 1 << std::endl;
      std::set< std::pair<uint32_t,float>> bad = BAD[APV];
      for(it=bad.begin();it!=bad.end();it++) (*oBAD[APV]) << (*it).first << "   " << 1  << std::endl;
      std::set< std::pair<uint32_t,float>> toBeChecked = TBC[APV];
      for(it=toBeChecked.begin();it!=toBeChecked.end();it++) (*oTBC[APV]) << (*it).first << "   " << 1 << std::endl;
      std::set< std::pair<uint32_t,float>> extreme = EXT[APV];
      for(it=extreme.begin();it!=extreme.end();it++) (*oEXT[APV]) << (*it).first << "   " << 1 << std::endl;
      std::set< std::pair<uint32_t,float>> noisy = NSY[APV];
      for(it=noisy.begin();it!=noisy.end();it++) (*oNSY[APV]) << (*it).first << "   " << 1 << std::endl;
    } 

    for (unsigned int i=0; i<3; i++) {
      oUNC[i]->close() ; delete oUNC[i];
      oBAD[i]->close() ; delete oBAD[i];
      oTBC[i]->close() ; delete oTBC[i];
      oEXT[i]->close() ; delete oEXT[i];
      oNSY[i]->close() ; delete oNSY[i];
    }
}


void SiStripApvGainFromFileBuilder::output_ratio_maps() const {
    // check if directory exist, if not create it
    struct stat sb;
    if ( !(stat("tracker_maps",&sb)==0 && S_ISDIR(sb.st_mode)) ) {
      if( mkdir("tracker_maps",0777)==-1 ) throw "cannot create the tracker_maps directory!";
    }

    // open output file
    std::ofstream* ofile_1 = new std::ofstream("tracker_maps/Ratio_APVpair1.txt", std::ofstream::trunc);
    if ( !ofile_1->is_open() ) throw "cannot open tracker_maps/Ratio_APVpair1.txt file!";
    std::ofstream* ofile_2 = new std::ofstream("tracker_maps/Ratio_APVpair2.txt", std::ofstream::trunc);
    if ( !ofile_2->is_open() ) throw "cannot open tracker_maps/Ratio_APVpair2.txt file!";
    std::ofstream* ofile_3 = new std::ofstream("tracker_maps/Ratio_APVpair3.txt", std::ofstream::trunc);
    if ( !ofile_3->is_open() ) throw "cannot open tracker_maps/Ratio_APVpair3.txt file!";


    for(unsigned int i=0; i<summary_.size(); i++) {
      Summary s = summary_[i];
      if( s.tickmark_height>0. && s.reference_height>0. && (s.aoh_gain-s.reference_aoh_gain)==0 ) {
        float ratio = s.tickmark_height / s.reference_height;
        if( s.offlineAPV_id==0 ) (*ofile_1) << s.det_id  << "    " << ratio << std::endl;
        else if ( s.offlineAPV_id==2 ) (*ofile_2) << s.det_id  << "    " << ratio << std::endl;
        else if ( s.offlineAPV_id==4 ) (*ofile_3) << s.det_id  << "    " << ratio << std::endl;
      }
    }

    for(unsigned int i=0; i<ex_summary_.size(); i++) {
      Summary s = ex_summary_[i];
      if( s.tickmark_height>0. && s.reference_height>0. && (s.aoh_gain-s.reference_aoh_gain)==0 ) {
        float ratio = s.tickmark_height / s.reference_height;
        if( s.offlineAPV_id==0 ) (*ofile_1) << s.det_id  << "    " << ratio << std::endl;
        else if ( s.offlineAPV_id==2 ) (*ofile_2) << s.det_id  << "    " << ratio << std::endl;
        else if ( s.offlineAPV_id==4 ) (*ofile_3) << s.det_id  << "    " << ratio << std::endl;
      }
    }

    ofile_1->close();
    delete ofile_1;

    ofile_2->close();
    delete ofile_2;

    ofile_3->close();
    delete ofile_3;
}

void SiStripApvGainFromFileBuilder::output_tickmark_maps(std::vector<Gain*>* maps, const char* basename) const {
    for (unsigned int APV=0; APV<maps->size(); APV+=2) {
        Gain* map = (*maps)[APV];
        if ( map!=0 ) {
            // check if directory exist, if not create it
            struct stat sb;
            if ( !(stat("tracker_maps",&sb)==0 && S_ISDIR(sb.st_mode)) ) {
              if( mkdir("tracker_maps",0777)==-1 ) throw "cannot create the tracker_maps directory!";
            }

            // open output file
            std::stringstream name;
            name << "tracker_maps/" << basename << "_APVpair" << (APV/2)+1 << ".txt";
            std::ofstream* ofile = new std::ofstream(name.str(), std::ofstream::trunc);
            if ( !ofile->is_open() ) throw "cannot open output file!";
            for(Gain::const_iterator el = map->begin(); el!=map->end(); el++) {
                (*ofile) << (*el).first << "    " << (*el).second.tickmark_height << std::endl;
            }
            ofile->close();
            delete ofile;
        }
    }
}


void SiStripApvGainFromFileBuilder::output_summary() const {

    std::string header = "  det id |APV| Connected |FED|CH|FEC|FECs|FECr|CCU|CCUc|i2cA|APVON|  noise  |tickHeig|AOH|Bias|DBgain|refHeigh|AOH|Bias|Status|";
    std::string ruler  = "---------+---+-----------+---+--+---+----+----+---+----+----+-----+---------+--------+---+----+------+--------+---+----+------+";

    // open output files
    std::ofstream* oSummary  = new std::ofstream("SiStripApvGainSummary.txt",std::ofstream::trunc);
    std::ofstream* oUncabled = new std::ofstream("SiStripApvGainUncabled.txt",std::ofstream::trunc);
    std::ofstream* oUncAndSc = new std::ofstream("SiStripApvGainUncAndSc.txt",std::ofstream::trunc);
    std::ofstream* oBad      = new std::ofstream("SiStripApvGainBad.txt",std::ofstream::trunc);
    std::ofstream* oNoisy    = new std::ofstream("SiStripApvGainNoisy.txt",std::ofstream::trunc);
    std::ofstream* oCheck    = new std::ofstream("SiStripApvGainToBeChecked.txt",std::ofstream::trunc);
    std::ofstream* oExtreme  = new std::ofstream("SiStripApvGainExtreme.txt",std::ofstream::trunc);

    (*oSummary)  << header << std::endl << ruler << std::endl;
    (*oUncabled) << header << std::endl << ruler << std::endl;
    (*oUncAndSc) << header << std::endl << ruler << std::endl;
    (*oBad)      << header << std::endl << ruler << std::endl;
    (*oNoisy)    << header << std::endl << ruler << std::endl;
    (*oCheck)    << header << std::endl << ruler << std::endl;
    (*oExtreme)  << header << std::endl << ruler << std::endl;


    for (unsigned int s=0; s<summary_.size(); s++) {
        Summary summary = summary_[s];
        std::stringstream line;
        format_summary(line, summary);

        if ( (!summary.is_connected) && (summary.is_scanned) ) (*oUncAndSc) << line.str() << std::endl;
        else if ( !summary.is_connected )                      (*oUncabled) << line.str() << std::endl;
        else if ( summary.tickmark_status=='B' )               (*oBad)      << line.str() << std::endl;
        else if ( summary.tickmark_status=='N' )               (*oNoisy)    << line.str() << std::endl;
        else if ( summary.tickmark_status=='V' )               (*oCheck)    << line.str() << std::endl;
        else if ( summary.tickmark_status=='E' )               (*oExtreme)  << line.str() << std::endl;
        else                                                   (*oSummary)  << line.str() << std::endl;      
    }

    for (unsigned int s=0; s<ex_summary_.size(); s++) {
        Summary summary = ex_summary_[s];
        std::stringstream line;
        format_summary(line, summary);

        if ( (!summary.is_connected) && (summary.is_scanned) ) (*oUncAndSc) << line.str() << std::endl;
        else if ( !summary.is_connected )                      (*oUncabled) << line.str() << std::endl;
        else if ( summary.tickmark_status=='B' )               (*oBad)      << line.str() << std::endl;
        else if ( summary.tickmark_status=='N' )               (*oNoisy)    << line.str() << std::endl;
        else if ( summary.tickmark_status=='V' )               (*oCheck)    << line.str() << std::endl;
        else if ( summary.tickmark_status=='E' )               (*oExtreme)  << line.str() << std::endl;
        else                                                   (*oSummary)  << line.str() << std::endl;
    }

    oSummary->close();   delete oSummary;
    oUncabled->close();  delete oUncabled;
    oUncAndSc->close();  delete oUncAndSc;
    oBad->close();       delete oBad;
    oNoisy->close();     delete oNoisy;
    oCheck->close();     delete oCheck;
    oExtreme->close();   delete oExtreme;
}


void SiStripApvGainFromFileBuilder::format_summary(std::stringstream& line, Summary summary) const {
    std::string conn = (summary.is_connected)? summary.detTopo()+" "+summary.detType()  : "NOT CONN";
    //std::string scan = (summary.is_scanned)?   "SCAN" : "NOT_SCAN";

    std::ios::fmtflags flags = line.flags();
    int tik_precision = (summary.tickmark_height<0 || summary.tickmark_height==999999.)? 0 : 3;
    int ref_precision = (summary.reference_height<0 || summary.reference_height==999999.)? 0 : 3;

    line << summary.det_id
         << std::setw(3)  << summary.offlineAPV_id << "  "
         << std::setw(11) << std::right << conn    << " " 
         << std::setw(3)  << summary.FED_id        << " " 
         << std::setw(2)  << summary.FED_ch        << "  " 
         << std::setw(2)  << summary.FEC_crate     << "  " 
         << std::setw(2)  << summary.FEC_slot      << "   " 
         << std::setw(2)  << summary.FEC_ring      << "  " 
         << std::setw(3)  << summary.CCU_addr      << "  " 
         << std::setw(2)  << summary.CCU_chan      << "  " 
         << std::setw(3)  << summary.i2cAdd        << "  " 
         << std::setw(3)  << summary.onlineAPV_id  << "   "
         << std::setw(8)  << summary.noise         << "   " 
         << std::setprecision(tik_precision) << std::fixed << std::setw(7) << summary.tickmark_height << " " 
         << std::setprecision(0) << std::setw(3) << summary.aoh_gain << "  " 
         << std::setw(3) << summary.aoh_bias
         << std::setprecision(3) << std::fixed << std::setw(7) << summary.gain_in_db << "  " 
         << std::setiosflags(flags) << std::setw(7) << std::setprecision(ref_precision)
         << summary.reference_height << " " 
         << std::setprecision(0) << std::setw(3) << summary.reference_aoh_gain << "  "
         << std::setw(3) << summary.reference_aoh_bias << " "
         << summary.tickmark_status  << " "
         << summary.channel_status << " "
         << summary.recovery_status << " "
         << summary.status_wrt_previous_run;
}


bool SiStripApvGainFromFileBuilder::height_from_maps(uint32_t det_id, int onlineAPV_id, float& tick_height) {

    Gain* map = 0;

    // search det_id and APV in the good scan map
    map = get_gain_map(&gains_, onlineAPV_id);
    if( map!=0 ) {
      Gain::const_iterator el = map->find( det_id );
      if ( el!=map->end() ) {
        tick_height = el->second.tickmark_height;
        return true;
      }
    }

    // search det_id and APV in the zero gain scan map
    map = get_gain_map(&negative_gains_, onlineAPV_id);
    if( map!=0 ) {
      Gain::const_iterator el = map->find( det_id );
      if ( el!=map->end() ) {
        tick_height = el->second.tickmark_height;
        return true;
      }
    }

    //search det_id and APV in the negative gain scan map
    map = get_gain_map(&null_gains_, onlineAPV_id);
    if( map!=0 ) {
      Gain::const_iterator el = map->find( det_id );
      if ( el!=map->end() ) {
        tick_height = el->second.tickmark_height;
        return true;
      }
    }

    return false; 
}

void SiStripApvGainFromFileBuilder::height_from_maps(uint32_t det_id, uint16_t totalAPVs, 
                                                   std::vector<std::pair<int,tickmark_info>>& tick_info) const {
  std::stringstream ex_msg;
  ex_msg << "two APVs with the same online id for det id " << det_id 
         << ". Please check the tick mark file or the read_tickmark routines." << std::endl;

  for (unsigned int i=0; i<6; i++) {
    int offlineAPV_id = online2offline(i, totalAPVs);
    try {
      Gain* map =  gains_.at(i);
      if( map!=0 ) {
        Gain::const_iterator el = map->find( det_id );
        if ( el!=map->end() ) {
	  if( tick_info[offlineAPV_id].second.tickmark_height!=999999. ) throw( ex_msg.str() );
          tick_info[offlineAPV_id].second=el->second;
        }
      }
    } catch (std::out_of_range) {
      // nothing to do, just pass over
    } 

    try {
      Gain* map =  negative_gains_.at(i);
      if( map!=0 ) {
        Gain::const_iterator el = map->find( det_id );
        if ( el!=map->end() ) {
          if( tick_info[offlineAPV_id].second.tickmark_height!=999999. ) throw( ex_msg.str() );
          tick_info[offlineAPV_id].second=el->second;
        }
      }
    } catch (std::out_of_range) {
      // nothing to do, just pass over
    }

    try {
      Gain* map =  null_gains_.at(i);
      if( map!=0 ) {
        Gain::const_iterator el = map->find( det_id );
        if ( el!=map->end() ) {
          if( tick_info[offlineAPV_id].second.tickmark_height!=999999. ) throw( ex_msg.str() );
          tick_info[offlineAPV_id].second=el->second;
        }
      }
    } catch (std::out_of_range) {
      // nothing to do, just pass over
    }
  }
}

void
SiStripApvGainFromFileBuilder::height_from_reference(uint32_t det_id, uint16_t totalAPVs, 
                                                     std::vector< std::pair<int,tickmark_info> >& tick_info) const {
  std::stringstream ex_msg;
  ex_msg << "two APVs with the same online id for det id " << det_id
         << ". Please check the tick mark file or the read_tickmark routines." << std::endl;

  for (unsigned int i=0; i<6; i++) {
    int offlineAPV_id = online2offline(i, totalAPVs);
    try {
      Gain* map =  ref_gains_.at(i);
      if( map!=0 ) {
        Gain::const_iterator el = map->find( det_id );
        if ( el!=map->end() ) {
          if( tick_info[offlineAPV_id].second.tickmark_height!=999999. ) throw( ex_msg.str() );
          tick_info[offlineAPV_id].second=el->second;
        }
      }
    } catch (std::out_of_range) {
      // nothing to do, just pass over
    }
  }
}

void
SiStripApvGainFromFileBuilder::noise_from_maps(uint32_t det_id, uint16_t totalAPVs,
                                               std::vector< std::pair<int,noise_info> >& noise_info) const {
  std::stringstream ex_msg;
  ex_msg << "two APVs with the same online id for det id " << det_id
         << ". Please check the noise file or the read_noise_from_root routine." << std::endl;

  for (unsigned int i=0; i<6; i++) {
    int offlineAPV_id = online2offline(i, totalAPVs);
    try {
      Noise* map =  noise_.at(i);
      if( map!=0 ) {
        Noise::const_iterator el = map->find( det_id );
        if ( el!=map->end() ) {
          if( noise_info[offlineAPV_id].second.noise!=999999. ) throw( ex_msg.str() );
          noise_info[offlineAPV_id].second=el->second;
        }
      }
    } catch (std::out_of_range) {
      // nothing to do, just pass over
    }
  }
}


void
SiStripApvGainFromFileBuilder::set_tickmark_status(Summary& summary) const {

  bool is_tickmark_usable  = is_usable( summary.tickmark_height );
  bool is_reference_usable = is_usable( summary.reference_height, summary.aoh_gain-summary.reference_aoh_gain);
  bool is_reference_ok     = is_usable( summary.reference_height );

  // Recovery procedure when the tickmark height is not bad or corrupted
  if ( is_tickmark_usable ) {

    if ( is_reference_usable ) {
      //summary.recovery_status = 'R';
      summary.status_wrt_previous_run = '=';
      float ratio = summary.tickmark_height/summary.reference_height;

      // check for channels with tickmark variation bigger than 5% w.r.t. reference
      if ( std::fabs(1.-ratio)>0.05 ) {
        summary.tickmark_status = 'V';
      }
    } else {
      //summary.recovery_status = 'U';
      if ( !is_reference_ok ) summary.status_wrt_previous_run = '+';
    }
  }

  // Recovery procedure for tickmark heights which are bad or corrupted
  if ( !is_tickmark_usable ) {
    // Recover if the receference height is usable
    if( is_reference_ok ) {
      summary.tickmark_status = 'B';
      //summary.recovery_status = 'R';
      summary.status_wrt_previous_run = '-';
    } else {
      summary.tickmark_status = 'B';
      //summary.recovery_status = 'U';
      summary.status_wrt_previous_run = '=';
    }
  }
}

float
SiStripApvGainFromFileBuilder::apply_recovery(Summary& summary) const {
  for(unsigned int i=0; i<recovery_.size(); i++) {
    if (summary.det_id==recovery_[i].det_id && summary.offlineAPV_id==recovery_[i].offlineAPV_id) {
      return recovery_[i].reference_height;
    }
  }
  return 0.;//summary.tickmark_height;
}

int SiStripApvGainFromFileBuilder::online2offline(uint16_t onlineAPV_id, uint16_t totalAPVs) const {
    return (onlineAPV_id>=totalAPVs)? onlineAPV_id -2 : onlineAPV_id;
}
