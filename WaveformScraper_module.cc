////////////////////////////////////////////////////////////////////////
// Class:       MCScraper
// Plugin Type: analyzer (Unknown Unknown)
// File:        MCScraper_module.cc
//
// Generated at Tue Jan 23 01:17:55 2024 by Alexander Antonakis using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "canvas/Persistency/Common/FindManyP.h"

// Additional LArSoft Includes
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "sbnobj/Common/Reco/OpT0FinderResult.h"
#include "sbnobj/Common/Reco/CRUMBSResult.h"
//#include "sbnobj/Common/Reco/CNNScore.h"
#include "sbnobj/Common/Reco/MVAPID.h"
#include "lardataobj/AnalysisBase/T0.h"

#include <cstring> // std::memcpy()
#include <vector>
#include <map>
#include <utility>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional> // std::mem_fn()
#include <typeinfo>
#include <cmath>

#include "TTimeStamp.h"

#include "art_root_io/TFileService.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <bitset>

#include "lardataobj/AnalysisBase/T0.h"
#include "larcoreobj/SummaryData/POTSummary.h"


// includes inspired from uboone code

#include "TDatabasePDG.h"
#include "TParticlePDG.h"



namespace analysis {
  class WaveformScraper;

  /// information from the subrun
  struct SubRunData_t {
    SubRunData_t() { Clear(); }
    void Clear() { pot = -99999.; }
    Double_t pot; //protons on target
    Double_t goodpot; // good pots
  }; // struct SubRunData_t

}


class analysis::WaveformScraper : public art::EDAnalyzer {
public:
  explicit WaveformScraper(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WaveformScraper(WaveformScraper const&) = delete;
  WaveformScraper(WaveformScraper&&) = delete;
  WaveformScraper& operator=(WaveformScraper const&) = delete;
  WaveformScraper& operator=(WaveformScraper&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginSubRun(const art::SubRun& sr);
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  //art::InputTag potSummaryTag_;
  // Declare member data here.

  // Create a SubRun TTree
  TTree *fSubRunTree;  

  SubRunData_t SubRunData;

  Double_t totPOT = 0.;
  Double_t totGoodPOT = 0.;
 
  TTree *fwaveTree; 


  // Tree variables
  unsigned int fEventID;
  unsigned int fRun;
  unsigned int fSubRun;

  double fEventTime;


  // Define waveform variables
  std::vector<raw::TimeStamp_t> fTimeStampVec;
  std::vector<raw::Channel_t> fChannelVec;
  //std::vector<std::vector<raw::ADC_Count_t>> fAdcVec;
  std::vector<std::vector<unsigned int>> fAdcVec;

  // Define input labels --> producers
  const std::string fPFParticleLabel; // Pandora for instance
  const std::string fTrackLabel;
  const std::string fShowerLabel; 
  const std::string fCalorimetryLabel; 
  const std::string fPOTModuleLabel;
  const std::string fVertexModuleLabel;
  const std::string fGenieGenLabel;
  const std::string fHitModuleLabel;
  const std::string fCosmicTaggerLabel;
  const std::string fOpT0FinderLabel;
  const std::string fCrumbsResultLabel;
  const std::string fPMTDecoderLabel;
  const std::string fPMTDecoderInstanceLabel;
 

  //const std::string fPMTDecoderLabel = "pmtdecoder";
  //const std::string fPMTDecoderInstanceLabel = "PMTChannels";
 //const std::string fPMTDecoderInstanceLabel = "FTrigChannels";

};


analysis::WaveformScraper::WaveformScraper(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fPFParticleLabel(p.get<std::string>("PFParticleLabel")), // Put all of the producer labels for DataProducts here
  fTrackLabel(p.get<std::string>("TrackLabel")),
  fShowerLabel(p.get<std::string>("ShowerLabel")),
  fCalorimetryLabel(p.get<std::string>("CalorimetryLabel")),
  fPOTModuleLabel(p.get<std::string>("POTModuleLabel")),
  fVertexModuleLabel(p.get<std::string>("VertexModuleLabel")),
  fGenieGenLabel(p.get<std::string>("GenieGenLabel")),
  fHitModuleLabel(p.get<std::string>("HitModuleLabel")),
  fCosmicTaggerLabel(p.get<std::string>("CosmicTaggerLabel")),
  fOpT0FinderLabel(p.get<std::string>("OpT0FinderLabel")),
  fCrumbsResultLabel(p.get<std::string>("CrumbsResultLabel")),
  fPMTDecoderLabel(p.get<std::string>("PMTDecoderLabel")),
  fPMTDecoderInstanceLabel(p.get<std::string>("PMTDecoderInstanceLabel"))
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}


void analysis::WaveformScraper::beginSubRun(const art::SubRun& sr)
{
  std::cout << "In the beginSubRun Function" << std::endl;
  art::Handle< sumdata::POTSummary > potListHandle;
  //sr.getByLabel(fPOTModuleLabel,potListHandle);

  if(sr.getByLabel(fPOTModuleLabel,potListHandle)) {
    SubRunData.pot=potListHandle->totpot;
    SubRunData.goodpot=potListHandle->totgoodpot;
    totPOT += potListHandle->totpot;
    totGoodPOT += potListHandle->totgoodpot;
    std::cout << "POT in SubRun " << potListHandle->totpot << std::endl;
    std::cout << "Good POT in SubRun " << potListHandle->totgoodpot << std::endl;
  }
  else {
    SubRunData.pot=0.;
  }
} // End of beginSubRun function



void analysis::WaveformScraper::analyze(art::Event const& e)
{
  // Grab the LArSoft Event ID for the current event
  fEventID = e.id().event();
  std::cout << "Processing Event: " << fEventID << std::endl;  

  fRun = e.run();
  fSubRun = e.subRun();
  std::cout << "Processing SubRun " << fSubRun << std::endl;

  fTimeStampVec.clear();
  fChannelVec.clear();
  fAdcVec.clear();

  art::Timestamp ts = e.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  fEventTime = tts.AsDouble();
  
  /*
  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  std::vector<art::Ptr<recob::PFParticle>> pfpVec;
  if (e.getByLabel(fPFParticleLabel, pfpHandle))
    art::fill_ptr_vector(pfpVec, pfpHandle);
  */

  

  // Grab the OpDetWaveforms
  art::Handle<std::vector<raw::OpDetWaveform>> waveHandle;
  std::vector<art::Ptr<raw::OpDetWaveform>> waveVec;
  if (e.getByLabel(fPMTDecoderLabel, fPMTDecoderInstanceLabel, waveHandle)) {
    std::cout << "pmtdecoder label is valid" << std::endl;
    art::fill_ptr_vector(waveVec, waveHandle);
  }
  else {
    std::cout << "PMT decoder label didn't work?" << std::endl;
  }
  // start loop over waveforms
  for (auto const &waveform : waveVec) {
  
    fTimeStampVec.push_back(waveform->TimeStamp());
    fChannelVec.push_back(waveform->ChannelNumber());
    //fAdcVec.push_back(waveform->Waveform());
    std::vector<raw::ADC_Count_t> curr_waveform = waveform->Waveform();
    //fAdcVec.push_back(curr_waveform);
    std::vector<unsigned int> temp_wave;
    std::cout << "Current waveform length " << curr_waveform.size() << std::endl;
    for (long unsigned int i = 0; i < curr_waveform.size(); ++i) {
      if (i < 5)
        std::cout << "w[i] = " << curr_waveform.at(i) << std::endl;
      temp_wave.push_back(curr_waveform.at(i));

    }
    fAdcVec.push_back(temp_wave);

  } // end loop over waveforms


  //copied from MergeDataPaddles.cxx
  //  art::Handle< raw::BeamInfo > beam;
  //    if (evt.getByLabel("beam",beam)){
  //        fData->beamtime = (double)beam->get_t_ms();
  //            fData->beamtime/=1000.; //in second
  //              }
  //

 
  //art::ServiceHandle<cheat::ParticleInventoryService> particleInventoryService;
  //auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);


  fwaveTree->Fill();

} // end of analyze event

void analysis::WaveformScraper::beginJob()
{
  
  // Implementation of optional member function here.
  std::cout << "/-------- Begin Job Step -------- /" << std::endl; 
  
  std::cout << std::endl;
  std::cout << "Clearing the SubRun Info at the Begin Job Step ???" << std::endl;
  SubRunData.Clear();

  art::ServiceHandle<art::TFileService> tfs;
 

  fSubRunTree = tfs->make<TTree>("subrun_ttree", "subrun_ttree");
  fSubRunTree->Branch("totPOT", &totPOT);
  fSubRunTree->Branch("totGoodPOT", &totGoodPOT);
 
  fwaveTree = tfs->make<TTree>("tree", "output ttree");
  fwaveTree->Branch("eventID", &fEventID);
  fwaveTree->Branch("Run", &fRun);
  fwaveTree->Branch("SubRun", &fSubRun);
  fwaveTree->Branch("EventTime", &fEventTime);
  fwaveTree->Branch("TimeStamp", &fTimeStampVec);
  fwaveTree->Branch("ChannelNumber", &fChannelVec);
  fwaveTree->Branch("ADC", &fAdcVec);

}

void analysis::WaveformScraper::endJob()
{
  // Implementation of optional member function here.

  // Fill the total POT Info at the end  
  fSubRunTree->Fill();
}

DEFINE_ART_MODULE(analysis::WaveformScraper)
