#ifndef SiStripApvGainFromFileBuilder_H
#define SiStripApvGainFromFileBuilder_H

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/ConditionDBWriter/interface/ConditionDBWriter.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "CondFormats/SiStripObjects/interface/SiStripApvGain.h"
#include "CalibFormats/SiStripObjects/interface/SiStripQuality.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#include <map>
#include <vector>
#include <stdint.h>


class SiStripDetCabling;
class RunInfo;
class TH1F;

class SiStripApvGainFromFileBuilder : public edm::EDAnalyzer {

 public:

  struct tickmark_info {
                   double tickmark_height;
                   double aoh_gain;
                   double aoh_bias;

                   tickmark_info(): tickmark_height(999999.), aoh_gain(256.), aoh_bias(256.) { }

                   tickmark_info(double height, double gain=256., double bias=256.) :
                                    tickmark_height(height), aoh_gain(gain), aoh_bias(bias) { }
                 };

  typedef std::map< uint32_t, tickmark_info > Gain;

  typedef struct {
                   uint32_t  det_id;
                   uint16_t  offlineAPV_id;
                   int       onlineAPV_id;
                   int       FED_id;
                   int       FED_ch;
                   int       i2cAdd;
                   int       FEC_crate;
                   int       FEC_slot;
                   int       FEC_ring;
                   int       CCU_addr;
                   int       CCU_chan;
                   bool      is_connected;
                   bool      is_scanned;
                   float     tickmark_height;
                   float     aoh_gain;
                   float     aoh_bias;
                   float     gain_in_db;
                   float     reference_height;
                   float     reference_aoh_gain;
                   float     reference_aoh_bias;
                   char      tickmark_status;           /*!< "G" good status; "B" bad status; "V" to be checked */
                   char      channel_status;            /*!< "e" enabled in FED; "d" disabled in FED */
                   char      recovery_status;           /*!< "R" recoverable; "U" unrecoverable (bad status in reference too) */
                   char      status_wrt_previous_run;   /*!< "+" better condition; "-" worse condition; "=" same condition */
                 } Summary;

  /** Brief Constructor.
   */ 
  explicit SiStripApvGainFromFileBuilder( const edm::ParameterSet& iConfig);

  /** Brief Destructor performing the memory cleanup.
   */ 
  ~SiStripApvGainFromFileBuilder();

  /** Brief One dummy-event analysis to create the database record.
   * This method loops over the full set of Det Id, inquires the internal 
   * tickmark maps to recover the associated heights - if any - plus other
   * infos (AOH gain and bias) and fill the gain into the DB payload. The
   * tickmark heights are translated into gains, dividing by 640. Exception
   * cases are identified and reported into the summary.
   */ 
  virtual void analyze(const edm::Event& , const edm::EventSetup& );

  /** Brief Routine to be called when the Job start.
   */ 
  virtual void beginJob();

  /** Brief Routine to be called when the Job end.
   */  
  virtual void endJob();


 private:
  edm::FileInPath gfp_;          /*!< File Path for the ideal geometry. */
  edm::FileInPath tfp_;          /*!< Path for the file containing the APV heights from the tickmark scan. */
  edm::FileInPath rfp_;          /*!< Path for the file containing the APV heights from the reference tickmark scan. */
  edm::FileInPath recoveryList_; /*!< Path for the file containing the list of channels to be recovered. */
  double heightThreshold_;       /*!< Lower threshold for accepting the APV tickmark height in the scan. */
  double goodHeightLimit_;       /*!< Lower threshold to consider the height values good for being used for recovery. */
  double badHeightLimit_;        /*!< Upper threshold to consider the height values bad and try the recovery procedure. */
  double dummyAPVGain_;          /*!< Dummy value for the APV gain. */
  bool doRecovery_;              /*!< Flag to apply the recovery procedure. */
  bool putDummyIntoUncabled_;    /*!< Flag for putting the dummy gain in the channels not actuall cabled. */
  bool putDummyIntoUnscanned_;   /*!< Flag for putting the dummy gain in the chennals not scanned. */
  bool putDummyIntoOffChannels_; /*!< Flag for putting the dummy gain in the channels that were off during the tickmark scan. */
  bool putDummyIntoBadChannels_; /*!< Flag for putting the dummy gain in the channels with negative tickmark heights. */
  bool putDummyIntoLowChannels_; /*!< Flag for putting the dummy gain in the channels with tickmark heights under threshold. */
  bool outputMaps_;              /*!< Flag for dumping the internal maps on ASCII files. */
  bool outputSummary_;           /*!< Flag for dumping the summary of the exceptions during the DB filling. */

  TH1F* h_tickmark_height_0;     /*!< Histogram storing the tickmark heights with AOH gain == 0. */
  TH1F* h_tickmark_height_1;     /*!< Histogram storing the tickmark heights with AOH gain == 1. */
  TH1F* h_tickmark_height_2;     /*!< Histogram storing the tickmark heights with AOH gain == 2. */
  TH1F* h_tickmark_height_3;     /*!< Histogram storing the tickmark heights with AOH gain == 3. */

  TH1F* h_tickmark_recovered_0;  /*!< Histogram storing the tickmark recovered heights with AOH gain == 0. */
  TH1F* h_tickmark_recovered_1;  /*!< Histogram storing the tickmark recovered heights with AOH gain == 1. */
  TH1F* h_tickmark_recovered_2;  /*!< Histogram storing the tickmark recovered heights with AOH gain == 2. */
  TH1F* h_tickmark_recovered_3;  /*!< Histogram storing the tickmark recovered heights with AOH gain == 3. */

  TH1F* h_tickmark_replaced_0;   /*!< Histogram storing the tickmark heights replaced with AOH gain == 0. */
  TH1F* h_tickmark_replaced_1;   /*!< Histogram storing the tickmark heights replaced with AOH gain == 1. */
  TH1F* h_tickmark_replaced_2;   /*!< Histogram storing the tickmark heights replaced with AOH gain == 2. */
  TH1F* h_tickmark_replaced_3;   /*!< Histogram storing the tickmark heights replaced with AOH gain == 3. */


  TH1F* h_height_ratio_at_DG0;   /*!< tickmark height / reference height for the channels with Delta AOH Gain ==  0. */
  TH1F* h_height_ratio_at_DGp1;  /*!< tickmark height / reference height for the channels with Delta AOH Gain == +1. */
  TH1F* h_height_ratio_at_DGp2;  /*!< tickmark height / reference height for the channels with Delta AOH Gain == +2. */
  TH1F* h_height_ratio_at_DGp3;  /*!< tickmark height / reference height for the channels with Delta AOH Gain == +3. */
  TH1F* h_height_ratio_at_DGm1;  /*!< tickmark height / reference height for the channels with Delta AOH Gain == -1. */
  TH1F* h_height_ratio_at_DGm2;  /*!< tickmark height / reference height for the channels with Delta AOH Gain == -2. */
  TH1F* h_height_ratio_at_DGm3;  /*!< tickmark height / reference height for the channels with Delta AOH Gain == -3. */


  edm::ESHandle<SiStripDetCabling> detCabling_;     /*!< Description of detector cabling. */
  SiStripQuality*                  siStripQuality_; /*!< Description of detector configuration (for the disabled channels). */

  /** Brief Maps [det_id <--> tickmark height infos] arranged per APV indexes.
   */
  std::vector<Gain*> gains_;          /*!< Mapping channels with positive heights. */
  std::vector<Gain*> negative_gains_; /*!< Mapping channels sending bad data. */
  std::vector<Gain*> null_gains_;     /*!< Mapping channels switched off during the scan. */

  std::vector<Gain*> ref_gains_;      /*!< Mapping channels with reference height for comparison. */


  /** Brief Collection of the channels entered in the DB without exceptions.
   * The channels whose APV gain has been input in the DB straight from the 
   * tickmark scan are collected in the summary vector. The summary list is
   * dumped in the SiStripApvGainSummary.txt at the end of the job. 
   */
  std::vector<Summary> summary_;    /*!< Collection of channel with no DB filling exceptions. */ 

  /** Brief Collection of the exceptions encountered when filling the DB. 
   * An exception occur for all the non-cabled channels ( no height associated
   * in the tikmark file) and for all the channels that were off ( with zero
   * height associated) or sending corrupted data (with negative values in the
   * tickmark file). At the end of the job the exception summary is dumped in
   * SiStripApvGainExceptionSummary.txt.
   */
  std::vector<Summary> ex_summary_; /*!< Collection of DB filling exceptions. */

  std::vector<Summary> recovery_;   /*!< Collection of channels to be recovered wiht the height from the reference run. */

  /** Brief Steers the reading of the summary file.
   * This method reads back a summary list and put it into the Summary
   * vector given in input.
   */
  bool read_summary(std::vector<Summary>& list);

  /** Brief Steers the reading of the tickmark heights.
   * This method clear the internal maps, cheks if the file is an ASCII file
   * or a root file and then calls the appropriate reading routine.
   */
  void read_tickmark(void);
  
  /** Brief Read the ASCII file containing the tickmark heights.
   * This method reads the ASCII files that contains the tickmark heights for 
   * every APV. Maps are created for every APV index. Negative and Zero heights,
   * yielding to a non physical gain, are stored into separate maps.
   *   Negative heights: channels sending bad data at the tickmark scan. 
   *   Zero heights    : channels switched off during the tickmark scan. 
   */ 
  bool read_tickmark_from_ascii(const char* filename, std::vector<Gain*>* gains, 
                                                      std::vector<Gain*>* negative_gains,
                                                      std::vector<Gain*>* null_gains);

  /** Brief Read the root file containing the tickmark heights.
   * This method reads the root file, generated by the WBM template, that
   * containis the tree with the tickmark height plus the values of the AOH
   * gain and bias. Maps with the tickmark heights are created for every APV
   * index. Negative and Zero heights, yielding to a non physical gain, are
   * stored into separate maps. 
   *   Negative heights: channels sending bad data at the tickmark scan. 
   *   Zero heights    : channels switched off during the tickmark scan. 
   */
  bool read_tickmark_from_root(const char* filename, std::vector<Gain*>* gains,
                                                     std::vector<Gain*>* negative_gains,
                                                     std::vector<Gain*>* null_gains);

  /** Brief Returns the mapping among det_id and tickmark heights for the APVs.
   * This method retrieves the mapping of detector Ids <-> tickmark heights for
   * a specific APV online Id. The map corresponding to the APV Id is returned 
   * if it exists; on the contrary a new empty map is created, inserted and 
   * returned. The methods accepts onlineIDs running from 0 to 5. 
   */
  Gain* get_map(std::vector<Gain*>* maps, int onlineAPV_id);

  /** Brief Dumps the tickmark heights in ASCII format to input the tracker map.
   * This method dumps the detector id <-> tickmark height into acii files for
   * each APV pair. The basename of for the acii file has to be provided as a
   * input parameter. These maps can be given to the tracker map tool
   */ 
  void output_tickmark_maps(std::vector<Gain*>* maps, const char* basename) const;

  /** Brief Dumps the the ratio of tickmark heights to input the tracker map.
   * This method dumps the detector id <-> ratios of tickmark heights of two
   * scan runs for each APV pair. The ratio points out big variation that can
   * be due to a hardware malfunctioning during the scan itself.
   */ 
  void output_ratio_maps() const;

  /** Brief Dumps the map of BAD, Unconnected and to be checked channels.
   */
  void output_topology_maps() const;

  /** Brief Dump the exceptions summary on a ASCII file.
   * This method dumps the online coordinate of the channels for which  there
   * was an exception for filling the database record. Exceptions are the non
   * cabled modules, the channels that were off during the tickmark scan, the
   * channels sending corrupted data duirng the tickmark scan. These exceptions
   * have been solved putting a dummy gain into the DB record or putting a zero
   * gain.
   */
  void output_summary() const;

  /** Brief Format the output line for the channel summary.
   */
  void format_summary (std::stringstream& line, Summary summary) const;

  /** Brief Find the tickmark height for a pair det_id, APV_id in the internal maps.
   */
  bool height_from_maps(uint32_t det_id, int onlineAPV_id, float& tick_height);
  void height_from_maps(uint32_t det_id, uint16_t totalAPVs, std::vector< std::pair<int,tickmark_info> >& tick_info) const;
  void height_from_reference(uint32_t det_id, uint16_t totalAPVs, std::vector< std::pair<int,tickmark_info> >& tick_info) const;

  /** Brief Check if tickmarck is recoverable.
   */
  bool is_recoverable(Summary& summary) const;  

  /** Brief Apply the recovery using the recovery list.
   * The method checks the input summary agains the recovery list. If a 
   * corresponding entry is found the reference height in the recovery
   * entry is returned. If not, the height of the input summary is returned.
   */
  float apply_recovery(Summary& summary) const; 

  /** Brief Convert online APV id into offline APV id.
   */
  int online2offline(uint16_t onlineAPV_id, uint16_t totalAPVs) const; 

  /** Brief Check if tickmark height is usable.
   */
  inline bool is_usable(double height, double deltaAOHGain=0.) const; 
};

inline bool
SiStripApvGainFromFileBuilder::is_usable(double height, double deltaAOHGain) const {
  return ( height!=999999. && height>heightThreshold_ && deltaAOHGain==0 );
}


#endif
