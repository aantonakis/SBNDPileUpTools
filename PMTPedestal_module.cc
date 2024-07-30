////////////////////////////////////////////////////////////////////////
// Class:       PMTPedestal
// Plugin Type: analyzer (Unknown Unknown)
// File:        PMTPedestal_module.cc
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
#include "TMath.h"


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

#include "helper_functions.h"


namespace analysis {
  class PMTPedestal;

}


class analysis::PMTPedestal : public art::EDAnalyzer {
public:
  explicit PMTPedestal(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PMTPedestal(PMTPedestal const&) = delete;
  PMTPedestal(PMTPedestal&&) = delete;
  PMTPedestal& operator=(PMTPedestal const&) = delete;
  PMTPedestal& operator=(PMTPedestal&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  //art::InputTag potSummaryTag_;
  // Declare member data here.
 
  unsigned int fEventID;

  TNtuple *fpedestalTree; 
  
  // store the channels and means as we loop through the events
  std::vector<int> channels;
  std::vector<int> means;

  // Define waveform variables
  std::vector<raw::TimeStamp_t> fTimeStampVec;
  std::vector<raw::Channel_t> fChannelVec;
  //std::vector<std::vector<raw::ADC_Count_t>> fAdcVec;
  std::vector<std::vector<unsigned int>> fAdcVec;

  // Define input labels --> producers
  
  const std::string fPMTDecoderLabel;
  const std::string fPMTDecoderInstanceLabel;
 

};


analysis::PMTPedestal::PMTPedestal(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fPMTDecoderLabel(p.get<std::string>("PMTDecoderLabel")),
  fPMTDecoderInstanceLabel(p.get<std::string>("PMTDecoderInstanceLabel"))
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}



void analysis::PMTPedestal::analyze(art::Event const& e)
{
  // Grab the LArSoft Event ID for the current event
  fEventID = e.id().event();
  std::cout << "Processing Event: " << fEventID << std::endl;  

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
  
    //fTimeStampVec.push_back(waveform->TimeStamp());
    //fChannelVec.push_back(waveform->ChannelNumber());

    unsigned int pmt_channel = waveform->ChannelNumber();
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

    // calculate the mean ADC of the current waveform
    int mean = static_cast<int>(TMath::Mean(temp_wave.begin(), temp_wave.end()));
    channels.push_back(pmt_channel);
    means.push_back(mean);
    //fAdcVec.push_back(temp_wave);

  } // end loop over waveforms

} // end of analyze event

void analysis::PMTPedestal::beginJob()
{
  
  // Implementation of optional member function here.
  std::cout << "/-------- Begin Job Step -------- /" << std::endl; 
  std::cout << std::endl;
 

  art::ServiceHandle<art::TFileService> tfs;
 
  fpedestalTree = tfs->make<TNtuple>("pedestal_tree", "output pmt pedestals", "Channel:Pedestal");
  //fwaveTree->Branch("eventID", &fEventID);
 

}

void analysis::PMTPedestal::endJob()
{
  // Implementation of optional member function here.

  // here we need to calculate the pedestals and fill our NTuple

  std::set<int> uniqueChannels(channels.begin(), channels.end());

  //std::map<int, int> channelPedestals;

  for (int elem : uniqueChannels) {
    std::vector<int> temp_means;
    std::vector<size_t> indices = findIndices(channels, elem);
    for (int ind : indices) {
      temp_means.push_back(means.at(ind));
    }
    double m = TMath::Mean(temp_means.begin(), temp_means.end());
    //channelPedestals[elem] = m;

    fpedestalTree->Fill(elem, m);

  }

}

DEFINE_ART_MODULE(analysis::PMTPedestal)
