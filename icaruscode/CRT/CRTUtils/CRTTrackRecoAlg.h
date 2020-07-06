#ifndef CRTTRACKRECOALG_H_SEEN
#define CRTTRACKRECOALG_H_SEEN

//////////////////////////////////////////////////////////////////////////////////
// CRTTrackRecoAlg.h
//
// Functions for CRT track reconstruction
// written by T Brooks (tbrooks@fnal.gov), November 2018
// ported to and modified for use with icaruscode by Chris.Hilgenberg@colostate.edu
///////////////////////////////////////////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
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
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()
#include "icaruscode/CRT/CRTProducts/CRTHit.hh"
#include "icaruscode/CRT/CRTProducts/CRTTrack.hh"
#include "icaruscode/CRT/CRTUtils/CRTHitRecoAlg.h"

// c++
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <cmath> 
#include <memory>

// ROOT
#include "TVector3.h"
#include "TGeoManager.h"

namespace icarus{

  class CRTTrackRecoAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> TimeLimit {
        Name("TimeLimit"),
        Comment("")
      };

      fhicl::Atom<double> AverageHitDistance {
        Name("AverageHitDistance"),
        Comment("Distance to average hits over on same plane")
      };

      fhicl::Atom<double> DistanceLimit {
        Name("DistanceLimit"),
        Comment("Distance to combine CRT hits into track")
      };

    };

    CRTTrackRecoAlg(const Config& config);

    CRTTrackRecoAlg(const fhicl::ParameterSet& pset) :
      CRTTrackRecoAlg(fhicl::Table<Config>(pset, {})()) {}

    CRTTrackRecoAlg(double aveHitDist, double distLim);

    ~CRTTrackRecoAlg();

    void reconfigure(const Config& config);

    std::vector<std::vector<art::Ptr<crt::CRTHit>>> CreateCRTTzeros(std::vector<art::Ptr<crt::CRTHit>>);

    // Function to make creating CRTTracks easier
    crt::CRTTrack FillCrtTrack(crt::CRTHit hit1, crt::CRTHit hit2, bool complete);

    // Function to average hits within a certain distance of each other
    std::vector<std::pair<crt::CRTHit, std::vector<int>>> AverageHits(std::vector<art::Ptr<crt::CRTHit>> hits, std::map<art::Ptr<crt::CRTHit>, int> hitIds);
    std::vector<crt::CRTHit> AverageHits(std::vector<art::Ptr<crt::CRTHit>> hits);

    // Take a list of hits and find average parameters
    crt::CRTHit DoAverage(std::vector<art::Ptr<crt::CRTHit>> hits);

    // Create CRTTracks from list of hits
    std::vector<std::pair<crt::CRTTrack, std::vector<int>>> CreateTracks(std::vector<std::pair<crt::CRTHit, std::vector<int>>> hits);
    std::vector<crt::CRTTrack> CreateTracks(std::vector<crt::CRTHit> hits);

    // Calculate the tagger crossing point of CRTTrack candidate
    TVector3 CrossPoint(crt::CRTHit hit, TVector3 start, TVector3 diff);

  private:

    geo::GeometryCore const* fGeometryService;

    double fTimeLimit;
    double fAverageHitDistance;
    double fDistanceLimit;

    crt::CRTHitRecoAlg hitAlg;

  };

}

#endif
