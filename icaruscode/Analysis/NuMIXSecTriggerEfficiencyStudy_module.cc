////////////////////////////////////////////////////////////////////////
// Class:       NuMIXSecTriggerEfficiencyStudy
// Plugin Type: analyzer (art v3_04_00)
// File:        NuMIXSecTriggerEfficiencyStudy_module.cc
//
// Generated at Sun Jan 26 22:13:22 2020 by Wesley Ketchum using cetskelgen
// from cetlib version v3_09_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
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
#include "larcorealg/GeoAlgo/GeoAlgo.h"
#include "lardataobj/RawData/TriggerData.h"

#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

//output to ntuple
#include "art_root_io/TFileService.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TVector3.h"

namespace ana {
  class NuMIXSecTriggerEfficiencyStudy;
}


class ana::NuMIXSecTriggerEfficiencyStudy : public art::EDAnalyzer {
public:
  explicit NuMIXSecTriggerEfficiencyStudy(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuMIXSecTriggerEfficiencyStudy(NuMIXSecTriggerEfficiencyStudy const&) = delete;
  NuMIXSecTriggerEfficiencyStudy(NuMIXSecTriggerEfficiencyStudy&&) = delete;
  NuMIXSecTriggerEfficiencyStudy& operator=(NuMIXSecTriggerEfficiencyStudy const&) = delete;
  NuMIXSecTriggerEfficiencyStudy& operator=(NuMIXSecTriggerEfficiencyStudy&&) = delete;

  void analyze(art::Event const& e) override;
  void beginJob() override;
  void endJob() override;

  bool IsSignal(const art::Ptr<simb::MCTruth>& mctruth,
                const std::vector<art::Ptr<simb::MCParticle> >& mcptls,
                std::vector<size_t>& matched_indices,
                float* treevar);
  bool IsInFV(double x, double y, double z);
  float ContainedLength(const TVector3 &v0, const TVector3 &v1,
                      const std::vector<geoalgo::AABox> &boxes);

private:

  TNtuple* fTrigEffTree;
  //TTree* fTrigEffTree;

  int PassTrigger;
  Float_t ENu;
  Float_t MuonP;
  Float_t LeadingProtonP;
  Float_t LeadingChargedPionP;

  std::vector<std::vector<geo::BoxBoundedGeo>> fTPCVolumes;
  std::vector<geo::BoxBoundedGeo> fActiveVolumes;

  std::string fMCTruthLabel;
  std::string fMCParticleLabel;
  std::string fTriggerLabel;

  bool fSelectSignal;
  bool fDoDebug;

  int NTotalEvent;
  int NSignalEvent;

};


ana::NuMIXSecTriggerEfficiencyStudy::NuMIXSecTriggerEfficiencyStudy(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{

  fMCTruthLabel = p.get<std::string>("MCTruthLabel", "generator");
  fMCParticleLabel = p.get<std::string>("MCParticleLabel", "largeant");
  fTriggerLabel = p.get<std::string>("TriggerLabel", "emuTrigger");
  fDoDebug = p.get<bool>("DoDebug", false);
  fSelectSignal = p.get<bool>("SelectSignal", true);

  NTotalEvent = 0;
  NSignalEvent = 0;

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

}

void ana::NuMIXSecTriggerEfficiencyStudy::beginJob()
{

  // get access to the TFile service
  art::ServiceHandle<art::TFileService> tfs;

/*
  fTrigEffTree = tfs->make<TTree>("TrigEffTree", "");

  fTrigEffTree->Branch("PassTrigger/I", &PassTrigger);
  fTrigEffTree->Branch("ENu/F", &ENu);
  fTrigEffTree->Branch("MuonP/F", &MuonP);
  fTrigEffTree->Branch("LeadingProtonP/F", &LeadingProtonP);
  fTrigEffTree->Branch("LeadingChargedPionP/F", &LeadingChargedPionP);
*/

  fTrigEffTree = tfs->make<TNtuple>(
    "TrigEffTree",
    "",
    "PassTrigger:ENu:MuonP:LeadingProtonP:LeadingChargedPionP:MuonLength:LeadingProtonLength:LeadingChargedPionLength"
  );

}

void ana::NuMIXSecTriggerEfficiencyStudy::endJob()
{

  printf("[NuMIXSecTriggerEfficiencyStudy::endJob] %d / %d\n", NSignalEvent, NTotalEvent);

  double NonSigFrac = (NTotalEvent-NSignalEvent)/NTotalEvent;

  printf("[NuMIXSecTriggerEfficiencyStudy::endJob] -> Missing %d (%1.2f%%) by signal selection\n", NTotalEvent-NSignalEvent, NonSigFrac);

}

void ana::NuMIXSecTriggerEfficiencyStudy::analyze(art::Event const& e)
{

  // count event
  NTotalEvent += 1;

  // Get handle
  art::Handle<std::vector<simb::MCTruth>> mcTruthHandle = e.getHandle<std::vector<simb::MCTruth>>(fMCTruthLabel);
  art::Handle<std::vector<simb::MCParticle>> mcParticleHandle = e.getHandle<std::vector<simb::MCParticle>>(fMCParticleLabel);
  art::Handle<std::vector<raw::Trigger>> trig_handle = e.getHandle<std::vector<raw::Trigger>>(fTriggerLabel);

  // fill vector

  std::vector<art::Ptr<simb::MCTruth> > mctruths;
  art::fill_ptr_vector(mctruths, mcTruthHandle);

  std::vector<art::Ptr<simb::MCParticle> > mcptls;
  art::fill_ptr_vector(mcptls, mcParticleHandle);

  // back tracker
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  // Check trigger
  if(fDoDebug) std::cout << "@@ Number of Trigger = " << trig_handle->size() << std::endl;

  auto global_trigger_det_time = trig_handle->at(0).TriggerTime();
  auto beam_gate_det_time = trig_handle->at(0).BeamGateTime();
  double diff_ts = global_trigger_det_time - beam_gate_det_time;
  if(fDoDebug) std::cout << "- diff_ts = " << diff_ts << std::endl;

  float treevar_PassTrigger = (diff_ts>-5) ? +1. : -1.;

  bool HasSignalEvent = false;

  if(fDoDebug) std::cout << "@@ Number of MCTruth = " << mctruths.size() << std::endl;
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

    // Get the matched G4-MCParticle

    std::vector<size_t> matched_mcptl_indices;
    for (size_t i_mc=0; i_mc<mcptls.size(); i_mc++){
      auto const& mcptl = mcptls.at(i_mc);
      std::string start_process = mcptl->Process();
      if(start_process!="primary") continue;

      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(mcptl->TrackId());
      if(truth.get()==mctruth.get()){
        matched_mcptl_indices.push_back(i_mc);
      }
    }


    const int NTreeVar = fTrigEffTree->GetNvar();
    float* tree_vars = new float[NTreeVar];
    tree_vars[0] = treevar_PassTrigger;
    tree_vars[1] = mctruth->GetNeutrino().Nu().EndMomentum().Energy();

/*
    PassTrigger = treevar_PassTrigger;
    ENu = mctruth->GetNeutrino().Nu().EndMomentum().Energy();
    MuonP = -2.;
    LeadingProtonP = -2.;
    LeadingChargedPionP = -2.;
*/

    bool issignal = IsSignal(mctruth, mcptls, matched_mcptl_indices, tree_vars);

    if(fDoDebug) std::cout << "--> issignal = " << issignal << std::endl;

    if(fSelectSignal && !issignal) continue;

    HasSignalEvent = true;

    if(fDoDebug) std::cout << "--> this event is filled" << std::endl;

    fTrigEffTree->Fill(tree_vars);

  }

  if(HasSignalEvent){
    NSignalEvent += 1;
  }

}

bool ana::NuMIXSecTriggerEfficiencyStudy::IsSignal(
                const art::Ptr<simb::MCTruth>& mctruth,
                const std::vector<art::Ptr<simb::MCParticle> >& mcptls,
                std::vector<size_t>& matched_indices,
                float* treevar){

  const simb::MCNeutrino& nu = mctruth->GetNeutrino();

  int NuPDG = nu.Nu().PdgCode();

  int IsCC = (!nu.CCNC()) && (nu.Mode() != simb::kWeakMix);

  const TVector3& NuPosition = nu.Nu().Position().Vect();
  bool isfv = IsInFV(NuPosition.X(), NuPosition.Y(), NuPosition.Z());

  unsigned int nMu(0), nP(0), nPi(0);
  unsigned int genieNPhotons(0), genieNMesons(0), genieNBaryonsAndPi0(0);

  double LeadingMuon_Momentum = -1.;
  double LeadingProton_Momentum = -1.;
  double LeadingChargedPion_Momentum = -1.;

  double LeadingMuon_TrueLength = -1.;
  double LeadingProton_TrueLength = -1.;
  double LeadingChargedPion_TrueLength = -1.;



  bool passMuonPCut = false;
  bool passProtonPCut = false;

  for(const auto& idx: matched_indices){
    const auto& mcptl = mcptls.at(idx);
    int this_ptl_pdg = mcptl->PdgCode();

    if(!mcptl->NumberTrajectoryPoints()) continue;

    double momentum = mcptl->Momentum().Vect().Mag();
    double genE = mcptl->Momentum().E();

    // Get the entry and exit points

    bool cont_tpc = mcptl->NumberTrajectoryPoints() > 0;
    bool contained = mcptl->NumberTrajectoryPoints() > 0;

    int entry_point = -1;
    int exit_point = -1;
    int cryostat_index = -1;
    int tpc_index = -1;
    for (unsigned j = 0; j < mcptl->NumberTrajectoryPoints(); j++) {
      for (unsigned i = 0; i < fActiveVolumes.size(); i++) {
        if (fActiveVolumes.at(i).ContainsPosition(mcptl->Position(j).Vect())) {
          entry_point = j;
          cryostat_index = i;
          break;
        }
      }
      if (entry_point != -1) break;
    }

    double g4truelength = 0.;

    // Use every trajectory point if possible
    if (entry_point >= 0) {

      // find tpc of the entry point
      std::vector<geo::BoxBoundedGeo> volumes = fTPCVolumes.at(cryostat_index);
      for (unsigned i = 0; i < volumes.size(); i++) {
        if (volumes[i].ContainsPosition(mcptl->Position(entry_point).Vect())) {
          tpc_index = i;
          cont_tpc = entry_point == 0;
          break;
        }
        contained = entry_point == 0;
      }
      if (tpc_index < 0) {
        cont_tpc = false;
      }

      // setup aa volumes too for length calc
      // Define the volume used for length calculation to be the cryostat volume in question
      std::vector<geoalgo::AABox> aa_volumes;
      const geo::BoxBoundedGeo &v = fActiveVolumes.at(cryostat_index);
      aa_volumes.emplace_back(v.MinX(), v.MinY(), v.MinZ(), v.MaxX(), v.MaxY(), v.MaxZ());

      // Get the length and determine if any point leaves the active volume

      // particle trajectory
      const simb::MCTrajectory &trajectory = mcptl->Trajectory();
      TVector3 pos = trajectory.Position(entry_point).Vect();

      bool crosses_tpc = false;

      for (unsigned i = entry_point+1; i < mcptl->NumberTrajectoryPoints(); i++) {
        TVector3 this_point = trajectory.Position(i).Vect();
        // get the exit point
        // update if particle is contained
        // check if particle has crossed TPC
        if (!crosses_tpc) {
          for (unsigned j = 0; j < volumes.size(); j++) {
            if (volumes[j].ContainsPosition(this_point) && tpc_index >= 0 && j != ((unsigned)tpc_index)) {
              crosses_tpc = true;
              break;
            }
          }
        }
        // check if particle has left tpc
        if (cont_tpc) {
          cont_tpc = volumes[tpc_index].ContainsPosition(this_point);
        }

        if (contained) {
          contained = fActiveVolumes.at(cryostat_index).ContainsPosition(this_point);
        }

        // update length
        g4truelength += ContainedLength(this_point, pos, aa_volumes);

        if (!fActiveVolumes.at(cryostat_index).ContainsPosition(this_point) && fActiveVolumes.at(cryostat_index).ContainsPosition(pos)) {
          exit_point = i-1;
        }

        pos = trajectory.Position(i).Vect();
      }
    }

    // if still exit_point is not found
    if (exit_point < 0 && entry_point >= 0) {
      exit_point = mcptl->NumberTrajectoryPoints() - 1;
    }

    if ( abs(this_ptl_pdg) == 13 ) {
      nMu+=1;
      if ( momentum > LeadingMuon_Momentum ) {
        LeadingMuon_Momentum = momentum;
        LeadingMuon_TrueLength = g4truelength;
        passMuonPCut = (momentum > 0.226);
      }
    }
    if ( abs(this_ptl_pdg) == 2212 ) {
      nP+=1;
      if ( momentum > LeadingProton_Momentum ) {
        LeadingProton_Momentum = momentum;
        LeadingProton_TrueLength = g4truelength;
        passProtonPCut = (momentum > 0.4 && momentum < 1.);
      }
    }
    if ( abs(this_ptl_pdg) == 211 ) {
      if ( momentum > LeadingChargedPion_Momentum ) {
        LeadingChargedPion_Momentum = momentum;
        LeadingChargedPion_TrueLength = g4truelength;
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



  } // END loop mcparticle indicies


  treevar[2] = LeadingMuon_Momentum;
  treevar[3] = LeadingProton_Momentum;
  treevar[4] = LeadingChargedPion_Momentum;

  treevar[5] = LeadingMuon_TrueLength;
  treevar[6] = LeadingProton_TrueLength;
  treevar[7] = LeadingChargedPion_TrueLength;

  if(fDoDebug){
    std::cout << "[IsSignal] @@ Called" << std::endl;
    std::cout << "[IsSignal] - NuPDG = " << NuPDG << std::endl;
    std::cout << "[IsSignal] - IsCC = " << IsCC << std::endl;
    std::cout << "[IsSignal] - IsInFV = " << isfv << std::endl;
    std::cout << "[IsSignal] - nMu = " << nMu << std::endl;
    std::cout << "[IsSignal] - nP = " << nP << std::endl;
    std::cout << "[IsSignal] - nPi = " << nPi << std::endl;
    std::cout << "[IsSignal] - genieNPhotons = " << genieNPhotons << std::endl;
    std::cout << "[IsSignal] - genieNMesons = " << genieNMesons << std::endl;
    std::cout << "[IsSignal] - genieNBaryonsAndPi0 = " << genieNBaryonsAndPi0 << std::endl;
    std::cout << "[IsSignal] - LeadingMuon_Momentum = " << LeadingMuon_Momentum << std::endl;
    std::cout << "[IsSignal] - LeadingProton_Momentum = " << LeadingProton_Momentum << std::endl;
  }

  if(abs(NuPDG)!=14) return false;
  if(!IsCC) return false;
  if(!isfv) return false;

  if ( nMu!=1 || nP==0 || nPi > 0 || genieNPhotons > 0 || genieNMesons > 0 || genieNBaryonsAndPi0 > 0 ) return false;
  if ( !passMuonPCut) return false;
  if ( !passProtonPCut ) return false;

  return true;

}

bool ana::NuMIXSecTriggerEfficiencyStudy::IsInFV(double x, double y, double z){

  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;

  return (( ( x < -61.94 - 25 && x > -358.49 + 25 ) ||
            ( x >  61.94 + 25 && x <  358.49 - 25 )) &&
          ( ( y > -181.86 + 25 && y < 134.96 - 25 ) &&
            ( z > -894.95 + 30 && z < 894.95 - 50 ) ));

}

float ana::NuMIXSecTriggerEfficiencyStudy::ContainedLength(const TVector3 &v0, const TVector3 &v1,
                       const std::vector<geoalgo::AABox> &boxes) {
  static const geoalgo::GeoAlgo algo;
  // if points are the same, return 0
  if ((v0 - v1).Mag() < 1e-6) return 0;

  // construct individual points
  geoalgo::Point_t p0(v0);
  geoalgo::Point_t p1(v1);

  // construct line segment
  geoalgo::LineSegment line(p0, p1);

  double length = 0;

  // total contained length is sum of lengths in all boxes
  // assuming they are non-overlapping
  for (auto const &box: boxes) {
    int n_contained = box.Contain(p0) + box.Contain(p1);
    // both points contained -- length is total length (also can break out of loop)
    if (n_contained == 2) {
      length = (v1 - v0).Mag();
      break;
    }
    // one contained -- have to find intersection point (which must exist)
    if (n_contained == 1) {
      auto intersections = algo.Intersection(line, box);
      // Because of floating point errors, it can sometimes happen
      // that there is 1 contained point but no "Intersections"
      // if one of the points is right on the edge
      if (intersections.size() == 0) {
        // determine which point is on the edge
        double tol = 1e-5;
        bool p0_edge = algo.SqDist(p0, box) < tol;
        bool p1_edge = algo.SqDist(p1, box) < tol;
        assert(p0_edge || p1_edge);
        // contained one is on edge -- can treat both as not contained
        //
        // In this case, no length
        if ((p0_edge && box.Contain(p0)) || (box.Contain(p1) && p1_edge))
          continue;
        // un-contaned one is on edge -- treat both as contained
        else if ((p0_edge && box.Contain(p1)) || (box.Contain(p0) && p1_edge)) {
    length = (v1 - v0).Mag();
    break;
        }
        else {
          assert(false); // bad
        }
      }
      // floating point errors can also falsely cause 2 intersection points
      //
      // in this case, one of the intersections must be very close to the
      // "contained" point, so the total contained length will be about
      // the same as the distance between the two intersection points
      else if (intersections.size() == 2) {
        length += (intersections.at(0).ToTLorentzVector().Vect() - intersections.at(1).ToTLorentzVector().Vect()).Mag();
        continue;
      }
      // "Correct"/ideal case -- 1 intersection point
      else if (intersections.size() == 1) {
        // get TVector at intersection point
        TVector3 int_tv(intersections.at(0).ToTLorentzVector().Vect());
        length += ( box.Contain(p0) ? (v0 - int_tv).Mag() : (v1 - int_tv).Mag() );
      }
      else assert(false); // bad
    }
    // none contained -- either must have zero or two intersections
    if (n_contained == 0) {
      auto intersections = algo.Intersection(line, box);
      if (!(intersections.size() == 0 || intersections.size() == 2)) {
        // more floating point error fixes...
        //
        // figure out which points are near the edge
        double tol = 1e-5;
        bool p0_edge = algo.SqDist(p0, box) < tol;
        bool p1_edge = algo.SqDist(p1, box) < tol;
        // and which points are near the intersection
        TVector3 vint = intersections.at(0).ToTLorentzVector().Vect();

        bool p0_int = (v0 - vint).Mag() < tol;
        bool p1_int = (v1 - vint).Mag() < tol;
        // exactly one of them should produce the intersection
        assert((p0_int && p0_edge) != (p1_int && p1_edge));
        // void variables when assert-ions are turned off
        (void) p0_int; (void) p1_int;

        // both close to edge -- full length is contained
        if (p0_edge && p1_edge) {
          length += (v0 - v1).Mag();
        }
        // otherwise -- one of them is not on an edge, no length is contained
        else {}
      }
      // assert(intersections.size() == 0 || intersections.size() == 2);
      else if (intersections.size() == 2) {
        TVector3 start(intersections.at(0).ToTLorentzVector().Vect());
        TVector3 end(intersections.at(1).ToTLorentzVector().Vect());
        length += (start - end).Mag();
      }
    }
  }

  return length;

}



DEFINE_ART_MODULE(ana::NuMIXSecTriggerEfficiencyStudy)
