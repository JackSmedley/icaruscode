////////////////////////////////////////////////////////////////////////
// Class:       NuMIXSecGenFilter
// Plugin Type: analyzer (art v3_04_00)
// File:        NuMIXSecGenFilter_module.cc
//
// Generated at Sun Jan 26 22:13:22 2020 by Wesley Ketchum using cetskelgen
// from cetlib version v3_09_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "lardataobj/RawData/TriggerData.h"

#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

//output to ntuple
#include "art_root_io/TFileService.h"
#include "TNtuple.h"

namespace ana {
  class NuMIXSecGenFilter;
}


class ana::NuMIXSecGenFilter : public art::EDFilter {
public:
  explicit NuMIXSecGenFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuMIXSecGenFilter(NuMIXSecGenFilter const&) = delete;
  NuMIXSecGenFilter(NuMIXSecGenFilter&&) = delete;
  NuMIXSecGenFilter& operator=(NuMIXSecGenFilter const&) = delete;
  NuMIXSecGenFilter& operator=(NuMIXSecGenFilter&&) = delete;

  virtual bool filter(art::Event& e) override;
  bool IsInFV(double x, double y, double z);
  void beginJob() override;

private:

/*
  std::vector<std::vector<geo::BoxBoundedGeo>> fTPCVolumes;
  std::vector<geo::BoxBoundedGeo> fActiveVolumes;
*/
  std::string fMCTruthLabel;
  int fSelectionMode;
  bool fDoDebug;

};


ana::NuMIXSecGenFilter::NuMIXSecGenFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}
{

  fMCTruthLabel = p.get<std::string>("MCTruthLabel", "generator");
  fSelectionMode = p.get<int>("SelectionMode");
  fDoDebug = p.get<bool>("DoDebug", false);

/*
  // geometry
  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();

  // first the TPC volumes
  for (auto const &cryo: geometry->Iterate<geo::CryostatGeo>()) {
    std::vector<geo::BoxBoundedGeo> this_tpc_volumes;
    for (auto const& TPC : geometry->Iterate<geo::TPCGeo>(cryo.ID())) {
      this_tpc_volumes.push_back(TPC.ActiveBoundingBox());
    }
     fTPCVolumes.push_back(std::move(this_tpc_volumes));
  }

  // then combine them into active volumes
  for (const std::vector<geo::BoxBoundedGeo> &tpcs: fTPCVolumes) {
    double XMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinX() < rhs.MinX(); })->MinX();
    double YMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinY() < rhs.MinY(); })->MinY();
    double ZMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinZ() < rhs.MinZ(); })->MinZ();

    double XMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxX() < rhs.MaxX(); })->MaxX();
    double YMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxY() < rhs.MaxY(); })->MaxY();
    double ZMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxZ() < rhs.MaxZ(); })->MaxZ();

    fActiveVolumes.emplace_back(XMin, XMax, YMin, YMax, ZMin, ZMax);
  }

*/
}

void ana::NuMIXSecGenFilter::beginJob()
{

}
bool ana::NuMIXSecGenFilter::filter(art::Event& e)
{

  std::cout << "[NuMIXSecGenFilter::filter] Called" << std::endl;

  // Get handle
  art::Handle<std::vector<simb::MCTruth>> mcTruthHandle = e.getHandle<std::vector<simb::MCTruth>>(fMCTruthLabel);

  // fill vector

  std::vector<art::Ptr<simb::MCTruth> > mctruths;
  art::fill_ptr_vector(mctruths, mcTruthHandle);

  if(fDoDebug) std::cout << "@@ Number of MCTruth = " << mctruths.size() << std::endl;

  bool PassSelection = false;

  for (size_t i_mc=0; i_mc<mctruths.size(); i_mc++) {

    auto const& mctruth = mctruths.at(i_mc);

    if(fDoDebug) std::cout << "- " << i_mc << "-th MCTruth" << std::endl;

    if (mctruth->NeutrinoSet()) {
      if(fDoDebug) std::cout << "  - Neuetrino set" << std::endl;
    }
    else{
      if(fDoDebug) std::cout << "  - Non-neutrino set (skipping)" << std::endl;
      continue;
    }

    const simb::MCNeutrino& nu = mctruth->GetNeutrino();
    bool InFV = IsInFV(nu.Nu().Vx(), nu.Nu().Vy(), nu.Nu().Vz());
    if(!InFV) continue;

    int NuPDG = nu.Nu().PdgCode();

    bool IsNC = nu.CCNC()  && (nu.Mode() != simb::kWeakMix);
    bool IsCC = (!nu.CCNC()) && (nu.Mode() != simb::kWeakMix);

    int NParticles = mctruth->NParticles();
    if(fDoDebug) std::cout << "  - Number of particles = " << NParticles << std::endl;


    unsigned int nMu(0), nP(0), nPi(0);
    unsigned int genieNPhotons(0), genieNMesons(0), genieNBaryonsAndPi0(0);
    double maxMomentumMuon = -1.;
    double maxMomentumP = -1.;
    double maxMomentumChargePion = -1.;
    bool passMuonPCut = false;
    bool passProtonPCut = false;

    for(int i_ptl=0; i_ptl<NParticles; i_ptl++){
      const simb::MCParticle& mcptl = mctruth->GetParticle(i_ptl);
      int this_ptl_pdg = mcptl.PdgCode();
      std::string start_process = mcptl.Process();
      int this_ptl_statuscode = mcptl.StatusCode();

      if(this_ptl_statuscode!=1) continue;

      double momentum = mcptl.Momentum().Vect().Mag();
      double genE = mcptl.Momentum().E();

      if(fDoDebug){
        printf("  - %d-th MCParticle\n", i_ptl);
        printf("    - pdg = %d\n", this_ptl_pdg);
        printf("    - start_process = %s\n", start_process.c_str());
        printf("    - Status = %d\n", this_ptl_statuscode);
        printf("    - P = %1.3f\n", momentum);
        printf("    - E = %1.3f\n", genE);
      }

      if ( abs(this_ptl_pdg) == 13 ) {
        nMu+=1;
        if ( momentum > maxMomentumMuon ) {
          maxMomentumMuon = momentum;
          passMuonPCut = (momentum > 0.226);
        }
      }
      if ( abs(this_ptl_pdg) == 2212 ) {
        nP+=1;
        if ( momentum > maxMomentumP ) {
          maxMomentumP = momentum;
          passProtonPCut = (momentum > 0.4 && momentum < 1.);
        }
      }
      if ( abs(this_ptl_pdg) == 211 ) {
        if ( momentum > maxMomentumChargePion ) {
          maxMomentumChargePion = momentum;
        }
      }

      if ( abs(this_ptl_pdg) == 111 || abs(this_ptl_pdg) == 211 ) nPi+=1;
      // CHECK A SIMILAR DEFINITION AS MINERVA FOR EXTRA REJECTION OF UNWANTED THINGS IN SIGNAL DEFN.
      if ( abs(this_ptl_pdg) == 22 && genE > 0.01 ) genieNPhotons+=1;
      else if ( abs(this_ptl_pdg) == 211 || abs(this_ptl_pdg) == 321 || abs(this_ptl_pdg) == 323 ||
                this_ptl_pdg == 111 || this_ptl_pdg == 130 || this_ptl_pdg == 310 || this_ptl_pdg == 311 ||
                this_ptl_pdg == 313 || abs(this_ptl_pdg) == 221 || abs(this_ptl_pdg) == 331 ) genieNMesons+=1;
      else if ( this_ptl_pdg == 3112 || this_ptl_pdg == 3122 || this_ptl_pdg == 3212 || this_ptl_pdg == 3222 ||
                this_ptl_pdg == 4112 || this_ptl_pdg == 4122 || this_ptl_pdg == 4212 || this_ptl_pdg == 4222 ||
                this_ptl_pdg == 411 || this_ptl_pdg == 421 || this_ptl_pdg == 111 ) genieNBaryonsAndPi0+=1;

    }

    bool Is1muNp0pi = 
          abs(NuPDG)==14 &&
          IsCC &&
          nMu==1 && nP>0 && nPi==0 && genieNPhotons==0 && genieNMesons==0 && genieNBaryonsAndPi0==0 &&
          passMuonPCut &&
          passProtonPCut;
    bool Is1muNpNpi =
          abs(NuPDG)==14 &&
          IsCC &&
          nMu==1 && nP>0 && nPi>0 &&
          passMuonPCut &&
          passProtonPCut;
    bool Is1muNpi =
          abs(NuPDG)==14 &&
          IsCC &&
          nMu==1 && nPi>0;
    bool IsNCEnergeticPion = 
          IsNC &&
          maxMomentumChargePion>0.190;

    if(fDoDebug){
      printf("  ====> Summary\n");
      if(IsCC) printf("  - CC\n");
      if(IsNC) printf("  - NC\n");
      printf("  - (nMu, nP, nPi) = (%d, %d, %d)\n", nMu, nP, nPi);
      printf("  - (nPhoron, nMesonm nBar) = (%d, %d, %d)\n", genieNPhotons, genieNMesons, genieNBaryonsAndPi0);
      printf("  - Muon Pmax = %1.3f\n", maxMomentumMuon);
      printf("  - Proton Pmax = %1.3f\n", maxMomentumP);
      printf("  - Charged pion Pmax = %1.3f\n", maxMomentumChargePion);
    }

    // CC 1muNp0pi
    if(fSelectionMode==0){
      if(Is1muNp0pi){
        PassSelection = true;
        if(fDoDebug) std::cout << "====> CC 1muNp0pi found" << std::endl;
      }
    }
    // Additional mode
    else if(fSelectionMode==2){
      if(Is1muNp0pi && maxMomentumMuon<0.4){
        PassSelection = true;
        if(fDoDebug) std::cout << "====> CC 1muNp0pi with muon P<0.4 found" << std::endl;
      }
    }

    // NC, at least one charged pion with ke>0.1GeV ~ p>0.190
    else if(fSelectionMode==1){
      if(IsNCEnergeticPion){
        PassSelection = true;
        if(fDoDebug) std::cout << "====> NC energetic pi found" << std::endl;
      }
    }

    // Is1muNpNpi
    else if(fSelectionMode==3){
      if(Is1muNpNpi){
        PassSelection = true;
        if(fDoDebug) std::cout << "====> Is1muNpNpi found" << std::endl;
      }
    }

    // Is1muNpi
    else if(fSelectionMode==4){
      if(Is1muNpi){
        PassSelection = true;
        if(fDoDebug) std::cout << "====> Is1muNpi found" << std::endl;
      }
    }

  }

  return PassSelection;

}

bool ana::NuMIXSecGenFilter::IsInFV(double x, double y, double z){

  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;

  return (( ( x < -61.94 - 25 && x > -358.49 + 25 ) ||
            ( x >  61.94 + 25 && x <  358.49 - 25 )) &&
          ( ( y > -181.86 + 25 && y < 134.96 - 25 ) &&
            ( z > -894.95 + 30 && z < 894.95 - 50 ) ));

}

DEFINE_ART_MODULE(ana::NuMIXSecGenFilter)
