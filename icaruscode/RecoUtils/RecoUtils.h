#ifndef FDSELECTIONUTILS_H_SEEN
#define FDSELECTIONUTILS_H_SEEN


///////////////////////////////////////////////
// RecoUtils.h
//
// A few reco utilities like truth matching 
// D Brailsford (adapted from work by D Brailsford and M Wallbank), October 2017
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h" 
#include "art_root_io/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
//#include "lardataobj/RecoBase/Track.h"
//#include "lardataobj/RecoBase/Shower.h"
//#include "lardataobj/AnalysisBase/MVAPIDResult.h"
//#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larcore/Geometry/Geometry.h"


// c++
#include <vector>
#include <map>

// ROOT
#include "TTree.h"


namespace RecoUtils{
  int TrueParticleID(const art::Ptr<recob::Hit>& hit); //Returns the geant4 ID which contributes the most to a single reco hit.  The matching method looks for true particle which deposits the most true energy in the reco hit
  int TrueParticleIDFromTotalTrueEnergy(const std::vector<art::Ptr<recob::Hit> >& hits); //Returns the geant4 ID which contributes the most to the vector of hits.  The matching method looks for which true particle deposits the most true energy in the reco hits
  int TrueParticleIDFromTotalRecoCharge(const std::vector<art::Ptr<recob::Hit> >& hits);  //Returns the geant4 ID which contributes the most to the vector of hits.  The matching method looks for which true particle contributes the most reconstructed charge to the hit selection (the reco charge of each hit is correlated with each maximally contributing true particle and summed)
  int TrueParticleIDFromTotalRecoHits(const std::vector<art::Ptr<recob::Hit> >& hits);  //Returns the geant4 ID which contributes the most to the vector of hits.  The matching method looks for which true particle maximally contributes to the most reco hits
  bool IsInsideTPC(TVector3 position, double distance_buffer); //Checks if a position is within any of the TPCs in the geometry (user can define some distance buffer from the TPC walls)
}

#endif
