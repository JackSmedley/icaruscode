////////////////////////////////////////////////////////////////////////
// Class:       TPCPMTMatchingAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        TPCPMTMatchingAna_module.cc
//
// Generated at Mon Jan 23 16:22:55 2023 by John Smedley using cetskelgen
// from  version .
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
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

//Data product includes
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"

//ROOT includes
#include "TTree.h"
#include "TVector3.h"



class TPCPMTMatchingAna;


class TPCPMTMatchingAna : public art::EDAnalyzer {
public:
  explicit TPCPMTMatchingAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TPCPMTMatchingAna(TPCPMTMatchingAna const&) = delete;
  TPCPMTMatchingAna(TPCPMTMatchingAna&&) = delete;
  TPCPMTMatchingAna& operator=(TPCPMTMatchingAna const&) = delete;
  TPCPMTMatchingAna& operator=(TPCPMTMatchingAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.
  void InitializeSlice();
  double CentroidOverlap(double center1, double center2, double width1, double width2);
  double CalculateAsymmetry(recob::OpFlash flash, int cryo);

  // Input parameters
  std::vector<std::string>  fCryoTags;             ///< Labels for each cryostat
  std::string               fOpFlashLabel;         ///< Label for PMT reconstruction products
  std::string               fPandoraLabel;         ///< Label for Pandora output products
  bool                      fCollectionOnly;       ///< Only use TPC spacepoints from the collection plane
  bool                      fVerbose;              ///< Print extra info, fcl input
  bool                      fFlashSelection;       ///< Apply Andrea Scarpeli's flash quality cuts
  int                       fMinimumPMTs;          ///< Used with fFlashSelection, number of PMTs with  >=1 OpHit over threshold required
  double                    fADCThreshold;         ///< Used with fFlashSelection, OpHit amplitude threshold (ADC)
  double                    fNominalTrigTime;      ///< Typical time of triggering flash, EYEBALLED (us)
  double                    fTriggerTolerance;     ///< Spread of triggering flash times, EYEBALLED (us)

   // Event-level data members
  int                       fRun;                  ///< Number of the run being processed
  int                       fEvent;                ///< Number of the event being processed
  int                       fCryo;                 ///< Cryostat this event occured in
  int                       fSliceNum;             ///< Number of slice in the event
  double                    fChargeT0;             ///< Start time for cathode-crossing PFPs, not always available (us)
  double                    fChargeTotal;          ///< Total charge in slice
  double                    fChargeCenterXGlobal;  ///< Weighted mean X position of spacepoints (cm)
  double                    fChargeCenterXLocal;   ///< Weighted mean X position of spacepoints, measured with respect to the cathode (cm)
  double                    fChargeCenterY;        ///< Weighted mean Y position of spacepoints (cm)
  double                    fChargeCenterZ;        ///< Weighted mean Z position of spacepoints (cm)
  double                    fChargeWidthX;         ///< Weighted standard deviation of X position of spacepoints (cm)
  double                    fChargeWidthY;         ///< Weighted standard deviation of Y position of spacepoints (cm)
  double                    fChargeWidthZ;         ///< Weighted standard deviation of Z position of spacepoints (cm)
  double                    fFlashT0;              ///< Matched OpFlash time OR earliest OpHit time associated to it, whichever value is less (us)
  double                    fFlashTime;            ///< Matched OpFlash time (us)
  double                    fFlashPEs;             ///< Brightness of matched flash (photoelectrons)
  double                    fFlashAsymmetry;       ///< East-West asymmetry of PEs in matched flash
  double                    fFlashCenterY;         ///< Weighted mean Y postion of hit PMTs (cm)
  double                    fFlashCenterZ;         ///< Weighted mean Z postion of hit PMTs (cm)
  double                    fFlashWidthY;          ///< Weighted standard deviation of Y postion of hit PMTs (cm)
  double                    fFlashWidthZ;          ///< Weighted standard deviation of Z postion of hit PMTs (cm)
  double                    fDeltaT;               ///< | Matched flash time - charge T0 | when available (us)
  double                    fDeltaY;               ///< | Matched flash Y center - charge Y center | (cm)
  double                    fDeltaZ;               ///< | Matched flash Z cetner - charge Z center | (cm)
  double                    fRadius;               ///< Hypotenuse of DeltaY and DeltaZ, PARAMETER MINIMIZED BY MATCHING (cm)
  double                    fOverlapY;             ///< Spacial overlap of flash and charge centroids in Y [>0] OR distance apart if no overlap [<0] (cm)
  double                    fOverlapZ;             ///< Spacial overlap of flash and charge centroids in Z [>0] OR distance apart if no overlap [<0] (cm)
  double                    fDeltaZ_Trigger;       ///< | Triggering flash Z cetner - charge Z center | (cm)

//Begin temporary TODO: Delete later

  double                    f2ndBestFlashT0;              ///< Matched OpFlash time OR earliest OpHit time associated to it, whichever value is less (us)
  double                    f2ndBestFlashTime;            ///< Matched OpFlash time (us)
  double                    f2ndBestFlashPEs;             ///< Brightness of matched flash (photoelectrons)
  double                    f2ndBestFlashAsymmetry;       ///< East-West asymmetry of PEs in matched flash
  double                    f2ndBestFlashCenterY;         ///< Weighted mean Y postion of hit PMTs (cm)
  double                    f2ndBestFlashCenterZ;         ///< Weighted mean Z postion of hit PMTs (cm)
  double                    f2ndBestFlashWidthY;          ///< Weighted standard deviation of Y postion of hit PMTs (cm)
  double                    f2ndBestFlashWidthZ;          ///< Weighted standard deviation of Z postion of hit PMTs (cm)
  double                    f2ndBestDeltaT;               ///< | Matched flash time - charge T0 | when available (us)
  double                    f2ndBestDeltaY;               ///< | Matched flash Y center - charge Y center | (cm)
  double                    f2ndBestDeltaZ;               ///< | Matched flash Z cetner - charge Z center | (cm)
  double                    f2ndBestRadius;               ///< Hypotenuse of DeltaY and DeltaZ, PARAMETER MINIMIZED BY MATCHING (cm)
  double                    f2ndBestOverlapY;             ///< Spacial overlap of flash and charge centroids in Y [>0] OR distance apart if no overlap [<0] (cm)
  double                    f2ndBestOverlapZ;             ///< Spacial overlap of flash and charge centroids in Z [>0] OR distance apart if no overlap [<0] (cm)

//End temporary

  TTree*                    fMatchTree;            ///< Tree to store all match information
};


TPCPMTMatchingAna::TPCPMTMatchingAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  // More initializers here.
  fCryoTags(p.get<std::vector<std::string>>("CryoTags")),
  fOpFlashLabel(p.get<std::string>("OpFlashLabel")),
  fPandoraLabel(p.get<std::string>("PandoraLabel")),
  fCollectionOnly(p.get<bool>("CollectionOnly")),
  fVerbose(p.get<bool>("Verbose")),
  fFlashSelection(p.get<bool>("FlashSelection")),
  fMinimumPMTs(p.get<int>("MinimumPMTs")),
  fADCThreshold(p.get<double>("ADCThreshold")),
  fNominalTrigTime(p.get<double>("NominalTrigTime")),
  fTriggerTolerance(p.get<double>("TriggerTolerance"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  art::ServiceHandle<art::TFileService> tfs;
  fMatchTree = tfs->make<TTree>("matchTree","TPC Slice - OpFlash Matching Analysis");

  //Event Info
  fMatchTree->Branch("run",                 &fRun,                 "run/I"                );
  fMatchTree->Branch("event",               &fEvent,               "event/I"              );
  fMatchTree->Branch("cryo",                &fCryo,                "cryo/I"               ); 
  fMatchTree->Branch("sliceNum",            &fSliceNum,            "sliceNum/I"           );

  //Charge Info
  fMatchTree->Branch("chargeT0",            &fChargeT0,            "chargeT0/d"           );
  fMatchTree->Branch("chargeTotal",         &fChargeTotal,         "chargeTotal/d"        );
  fMatchTree->Branch("chargeCenterXGlobal", &fChargeCenterXGlobal, "chargeCenterXGlobal/d");
  fMatchTree->Branch("chargeCenterXLocal",  &fChargeCenterXLocal,  "chargeCenterXLocal/d" );
  fMatchTree->Branch("chargeCenterY",       &fChargeCenterY,       "chargeCenterY/d"      );
  fMatchTree->Branch("chargeCenterZ",       &fChargeCenterZ,       "chargeCenterZ/d"      );
  fMatchTree->Branch("chargeWidthX",        &fChargeWidthX,        "chargeWidthX/d"       );
  fMatchTree->Branch("chargeWidthY",        &fChargeWidthY,        "chargeWidthY/d"       );
  fMatchTree->Branch("chargeWidthZ",        &fChargeWidthZ,        "chargeWidthZ/d"       );

  //Matched Flash Info
  fMatchTree->Branch("flashT0",             &fFlashT0,             "flashT0/d"            );
  fMatchTree->Branch("flashTime",           &fFlashTime,           "flashTime/d"          );
  fMatchTree->Branch("flashPEs",            &fFlashPEs,            "flashPEs/d"           );
  fMatchTree->Branch("flashAsymmetry",      &fFlashAsymmetry,      "flashAsymmetry/d"     );
  fMatchTree->Branch("flashCenterY",        &fFlashCenterY,        "flashCenterY/d"       );
  fMatchTree->Branch("flashCenterZ",        &fFlashCenterZ,        "flashCenterZ/d"       );
  fMatchTree->Branch("flashWidthY",         &fFlashWidthY,         "flashWidthY/d"        );
  fMatchTree->Branch("flashWidthZ",         &fFlashWidthZ,         "flashWidthZ/d"        );

  //Match Quality Info
  fMatchTree->Branch("deltaT",              &fDeltaT,              "deltaT/d"             );
  fMatchTree->Branch("deltaY",              &fDeltaY,              "deltaY/d"             );
  fMatchTree->Branch("deltaZ",              &fDeltaZ,              "deltaZ/d"             );
  fMatchTree->Branch("radius",              &fRadius,              "radius/d"             );
  fMatchTree->Branch("overlapY",            &fOverlapY,            "overlapY/d"           );
  fMatchTree->Branch("overlapZ",            &fOverlapZ,            "overlapZ/d"           );
  fMatchTree->Branch("deltaZ_Trigger",      &fDeltaZ_Trigger,      "deltaZ_Trigger/d"     );

//Begin temporary TODO: Delete later
  fMatchTree->Branch("2ndBestFlashT0",             &f2ndBestFlashT0,             "2ndBestFlashT0/d"            );
  fMatchTree->Branch("2ndBestFlashTime",           &f2ndBestFlashTime,           "2ndBestFlashTime/d"          );
  fMatchTree->Branch("2ndBestFlashPEs",            &f2ndBestFlashPEs,            "2ndBestFlashPEs/d"           );
  fMatchTree->Branch("2ndBestFlashAsymmetry",      &f2ndBestFlashAsymmetry,      "2ndBestFlashAsymmetry/d"     );
  fMatchTree->Branch("2ndBestFlashCenterY",        &f2ndBestFlashCenterY,        "2ndBestFlashCenterY/d"       );
  fMatchTree->Branch("2ndBestFlashCenterZ",        &f2ndBestFlashCenterZ,        "2ndBestFlashCenterZ/d"       );
  fMatchTree->Branch("2ndBestFlashWidthY",         &f2ndBestFlashWidthY,         "2ndBestFlashWidthY/d"        );
  fMatchTree->Branch("2ndBestFlashWidthZ",         &f2ndBestFlashWidthZ,         "2ndBestFlashWidthZ/d"        );

  fMatchTree->Branch("2ndBestDeltaT",              &f2ndBestDeltaT,              "2ndBestDeltaT/d"             );
  fMatchTree->Branch("2ndBestDeltaY",              &f2ndBestDeltaY,              "2ndBestDeltaY/d"             );
  fMatchTree->Branch("2ndBestDeltaZ",              &f2ndBestDeltaZ,              "2ndBestDeltaZ/d"             );
  fMatchTree->Branch("2ndBestRadius",              &f2ndBestRadius,              "2ndBestRadius/d"             );
  fMatchTree->Branch("2ndBestOverlapY",            &f2ndBestOverlapY,            "2ndBestOverlapY/d"           );
  fMatchTree->Branch("2ndBestOverlapZ",            &f2ndBestOverlapZ,            "2ndBestOverlapZ/d"           );
//End temporary

}

void TPCPMTMatchingAna::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  fEvent  = e.id().event();
  fRun    = e.run();

  //For cryo...
  for ( const std::string cryoTag : fCryoTags ) {
    //East-->0, West-->1
    fCryo = ( cryoTag.find("W") != std::string::npos ) ? 1 : 0;


/* ~~~~~~~~~~~~~~~~~~~~ Flash Section
 * Here we gather the OpFlashes found in this cryostat and their OpHits
 * We iterate through the flashes to identify a triggering flash and the earliest time associated with each flash
 * Additionally, we optionally apply Andrea's flash selection used for CRT-PMT studies (at least 5 hits in flash with amplitude >= 400)
 */

    //Fetch the flashes and their associated hits
    art::Handle<std::vector<recob::OpFlash>> flashHandle;
    e.getByLabel(fOpFlashLabel + cryoTag, flashHandle);
    art::FindMany<recob::OpHit> fmOpHits(flashHandle, e, fOpFlashLabel + cryoTag);

    int nFlashes = (*flashHandle).size();
    double triggerFlashCenter = -9999.;

    int nValidHits;
    double thisFlashTime;
    double minTime;
    std::vector<double> firstTimes;
    std::vector<bool> validFlashes;

    //For flash...
    for ( int i = 0; i < nFlashes; i++ ) {
      const recob::OpFlash &flash = (*flashHandle).at(i);
      const std::vector<recob::OpHit const*> &opHitsVec = fmOpHits.at(i);

      nValidHits = 0;
      thisFlashTime = flash.Time();
      minTime = 1e6;

      //Is this a triggering flash?
      if ( (thisFlashTime < fNominalTrigTime + fTriggerTolerance) && (thisFlashTime > fNominalTrigTime - fTriggerTolerance) ) triggerFlashCenter = flash.ZCenter();

      //For OpHit...
      for (const recob::OpHit *opHit : opHitsVec ) {
        if ( opHit->PeakTime() < minTime ) minTime = opHit->PeakTime();
        if ( fFlashSelection && opHit->Amplitude() >= fADCThreshold ) nValidHits++;
      } //End for OpHit

      firstTimes.push_back( std::min(minTime, thisFlashTime) ); //TODO: Understand why OpFlash.Time() is sometimes >2ns earlier than earliest hit
      validFlashes.push_back(nValidHits >= fMinimumPMTs); //TODO: Check number of PMTs with valid hits, not number of valid hits. "More correct" but maybe the same

      if ( minTime > thisFlashTime ) std::cout << "Cryo index " << cryoTag << " flash time " << thisFlashTime << " first hit time " << minTime << std::endl; //TODO: Remove print statements later

    } // End for flash

    std::cout << "Event: " << fEvent << ", Cryo: " << cryoTag << ", nFlashes: " << nFlashes << std::endl; //TODO: Remove print statements later


/* ~~~~~~~~~~~~~~~~~~~~ TPC Section
 * Here we start by gathering the Slices in the event
 * For each slice, the charge centroid is first calculated
 * Then we iterate through flashes to identify the best match flash
 * If a triggering flash was found earlier, the barycetner distance to the triggering flash is also stored
 */

    //Fetch slices and TPC hits
    art::Handle<std::vector<recob::Slice>> sliceHandle;
    e.getByLabel(fPandoraLabel + cryoTag, sliceHandle);
    art::FindManyP<recob::Hit> fmTPCHits(sliceHandle, e, fPandoraLabel + cryoTag);
    art::FindManyP<recob::PFParticle> fmPFPs(sliceHandle, e, fPandoraLabel + cryoTag);

    int nSlices = (*sliceHandle).size();

    //For slice...
    for ( int j = 0; j < nSlices; j++ ) {
      fSliceNum = j;
      InitializeSlice();

      const std::vector<art::Ptr<recob::Hit>> &tpcHitsVec = fmTPCHits.at(j);
      const std::vector<art::Ptr<recob::PFParticle>> &pfpsVec = fmPFPs.at(j);
      art::FindOne<recob::SpacePoint> f1SpacePoint(tpcHitsVec, e, fPandoraLabel + cryoTag);

      int nHits = tpcHitsVec.size();
      int nPFPs = pfpsVec.size();

      //These slices don't make it into the CAFs anyway, so why bother?
      if ( nPFPs == 1 || nPFPs == 0) {
        fMatchTree->Fill();
        if ( nPFPs == 1 ) std::cout << "1 PFP found in Event " << fEvent << " Slice " << j << ", PDG Code: " << pfpsVec.at(0)->PdgCode() << ", Continuing..." << std::endl; //TODO: Remove print statements later
        else std::cout << "0 PFPs found in Event " << fEvent << " Slice " << j << ", Continuing..." << std::endl; //TODO: Remove print statements later
        continue;
      }

      //Retrieve Pandora's T0 for this slice if available
      //Same for every PFP in slice so we only need one
      art::FindOne<anab::T0> f1T0( {pfpsVec.at(0)}, e, fPandoraLabel + cryoTag);
      if ( f1T0.at(0).isValid() ) fChargeT0 = f1T0.at(0).ref().Time() / 1e3;

      double sumCharge = 0.;
      double sumX = 0.; double sumY = 0.; double sumZ = 0.;
      double sumXX = 0.;double sumYY = 0.; double sumZZ = 0.;
      double thisHitCharge, thisHitXYZ[3];

      //For hit...
      for ( int k = 0; k < nHits; k++ ) {
        const art::Ptr<recob::Hit> &tpcHit = tpcHitsVec.at(k);

        //Only use hits with associated SpacePoints, and optionally only collection plane hits
        if ( fCollectionOnly && tpcHit->SignalType() != geo::kCollection ) continue;
        if ( !f1SpacePoint.at(k).isValid() ) continue; //TODO: Validate that we aren't missing boatloads of charge
        const recob::SpacePoint point = f1SpacePoint.at(k).ref();

        thisHitCharge = tpcHit->Integral();
        thisHitXYZ[0] = point.XYZ()[0];
        thisHitXYZ[1] = point.XYZ()[1];
        thisHitXYZ[2] = point.XYZ()[2];

        sumCharge += thisHitCharge;
        sumX += thisHitXYZ[0] * thisHitCharge;
        sumY += thisHitXYZ[1] * thisHitCharge;
        sumZ += thisHitXYZ[2] * thisHitCharge;
        sumXX += thisHitXYZ[0] * thisHitXYZ[0] * thisHitCharge;
        sumYY += thisHitXYZ[1] * thisHitXYZ[1] * thisHitCharge;
        sumZZ += thisHitXYZ[2] * thisHitXYZ[2] * thisHitCharge;
      } //End for hit

      //If all went correctly, this should never happen...
      if ( sumCharge == 0. ) {
        fMatchTree->Fill();
        std::cout << "No charge found in Event " << fEvent << " Slice " << j << "! Continuing..." << std::endl; //TODO: Remove print statements later, or maybe keep this one?
        continue;
      }

      //Update charge variables
      fChargeTotal = sumCharge;
      fChargeCenterXGlobal = sumX / sumCharge;
      //TODO: Get the cathode position and shift global X to local X in a less hacky way
      //According to a geometrydump, the cathode X positions are +/-(210.14, 210.29), depending on the TPC. Here I just averaged those...
      fChargeCenterXLocal = fChargeCenterXGlobal - 210.215 * (2*fCryo - 1);
      fChargeCenterY = sumY / sumCharge;
      fChargeCenterZ = sumZ / sumCharge;
      fChargeWidthX = pow( (sumXX/sumCharge - pow(sumX/sumCharge, 2)), .5);
      fChargeWidthY = pow( (sumYY/sumCharge - pow(sumY/sumCharge, 2)), .5);
      fChargeWidthZ = pow( (sumZZ/sumCharge - pow(sumZ/sumCharge, 2)), .5); 
      if ( triggerFlashCenter != -9999 ) fDeltaZ_Trigger = abs(triggerFlashCenter - fChargeCenterZ); 

      int matchIndex = -5;
      int scndBestMatchIndex = -5; //TODO:Delete later
      double minDistance = 1e6;
      double scndBestMinDistance = 1e6; //TODO:Delete later
      double thisFlashCenterY, thisFlashCenterZ, thisDistance;

      //For flash...
      for ( int m = 0; m < nFlashes; m++ ) {

        //Skip over flashes that didn't pass selection
        if ( fFlashSelection && validFlashes.at(m) == false ) continue;

        const recob::OpFlash &flash = (*flashHandle).at(m);
        thisFlashCenterY = flash.YCenter();
        thisFlashCenterZ = flash.ZCenter();

        //Find index of flash that minimizes barycenter distance in YZ place
        thisDistance = pow( (pow(thisFlashCenterY - fChargeCenterY, 2) + pow(thisFlashCenterZ - fChargeCenterZ, 2)), .5);
        if ( thisDistance < minDistance ) {
          scndBestMinDistance = minDistance; //TODO: Delete later
          scndBestMatchIndex = matchIndex; //TODO: Delete later
          minDistance = thisDistance;
          matchIndex = m;
        }
        else if ( thisDistance < scndBestMinDistance ) { //TODO: Delete later
          scndBestMinDistance = thisDistance;
          scndBestMatchIndex = m;
        }
      } //End for flash

      //If all went correctly, this should never happen...
      if ( matchIndex == -5 ) {
        fMatchTree->Fill();
        std::cout << "No matching flash found for Event " << fEvent << " Slice " << j << "! Continuing..." << std::endl; //TODO: Remove print statements later, or maybe keep this one?
        continue;
      }

      //Update match variables
      const recob::OpFlash &matchedFlash = (*flashHandle).at(matchIndex);
      double matchedTime = matchedFlash.Time();
      double matchedYCenter = matchedFlash.YCenter();
      double matchedZCenter = matchedFlash.ZCenter();
      double matchedYWidth = matchedFlash.YWidth();
      double matchedZWidth = matchedFlash.ZWidth();

      fFlashT0 = firstTimes.at(matchIndex);
      fFlashTime = matchedTime;
      fFlashPEs =  matchedFlash.TotalPE();
      fFlashAsymmetry = CalculateAsymmetry(matchedFlash, fCryo);
      fFlashCenterY = matchedYCenter;
      fFlashCenterZ = matchedZCenter;
      fFlashWidthY = matchedYWidth;
      fFlashWidthZ = matchedZWidth;
      if ( fChargeT0 !=-9999 ) fDeltaT = abs(matchedTime - fChargeT0);
      fDeltaY = abs(matchedYCenter - fChargeCenterY);
      fDeltaZ = abs(matchedZCenter - fChargeCenterZ);
      fRadius = pow( (pow(fDeltaY, 2) + pow(fDeltaZ, 2)), .5);
      fOverlapY = CentroidOverlap(matchedYCenter, fChargeCenterY, matchedYWidth, fChargeWidthY);
      fOverlapZ = CentroidOverlap(matchedZCenter, fChargeCenterZ, matchedZWidth, fChargeWidthZ);

      //Begin temporary TODO: Delete later
      if ( scndBestMatchIndex != -5 ) {

      const recob::OpFlash &scndBestMatchedFlash = (*flashHandle).at(scndBestMatchIndex);
      double scndBestMatchedTime = scndBestMatchedFlash.Time();
      double scndBestMatchedYCenter = scndBestMatchedFlash.YCenter();
      double scndBestMatchedZCenter = scndBestMatchedFlash.ZCenter();
      double scndBestMatchedYWidth = scndBestMatchedFlash.YWidth();
      double scndBestMatchedZWidth = scndBestMatchedFlash.ZWidth();

      f2ndBestFlashT0 = firstTimes.at(scndBestMatchIndex);
      f2ndBestFlashTime = scndBestMatchedTime;
      f2ndBestFlashPEs =  scndBestMatchedFlash.TotalPE();
      f2ndBestFlashAsymmetry = CalculateAsymmetry(scndBestMatchedFlash, fCryo);
      f2ndBestFlashCenterY = scndBestMatchedYCenter;
      f2ndBestFlashCenterZ = scndBestMatchedZCenter;
      f2ndBestFlashWidthY = scndBestMatchedYWidth;
      f2ndBestFlashWidthZ = scndBestMatchedZWidth;
      if ( fChargeT0 !=-9999 ) f2ndBestDeltaT = abs(scndBestMatchedTime - fChargeT0);
      f2ndBestDeltaY = abs(scndBestMatchedYCenter - fChargeCenterY);
      f2ndBestDeltaZ = abs(scndBestMatchedZCenter - fChargeCenterZ);
      f2ndBestRadius = pow( (pow(f2ndBestDeltaY, 2) + pow(f2ndBestDeltaZ, 2)), .5);
      f2ndBestOverlapY = CentroidOverlap(scndBestMatchedYCenter, fChargeCenterY, scndBestMatchedYWidth, fChargeWidthY);
      f2ndBestOverlapZ = CentroidOverlap(scndBestMatchedZCenter, fChargeCenterZ, scndBestMatchedZWidth, fChargeWidthZ);

      }
      //End temporary

      fMatchTree->Fill();

    } //End for slice

  } //End for cryo


} //End analyze()


void TPCPMTMatchingAna::InitializeSlice() {
  fChargeT0 = -9999.;
  fChargeTotal = -9999.;
  fChargeCenterXGlobal = -9999.;
  fChargeCenterXLocal = -9999.;
  fChargeCenterY = -9999.;
  fChargeCenterZ = -9999.;
  fChargeWidthX = -9999.;
  fChargeWidthY = -9999.;
  fChargeWidthZ = -9999.;
  fFlashT0 = -9999.;
  fFlashTime = -9999.;
  fFlashPEs = -9999.;
  fFlashCenterY = -9999.;
  fFlashCenterZ = -9999.;
  fFlashWidthY = -9999.;
  fFlashWidthZ = -9999.;
  fDeltaT = -9999.;
  fDeltaY = -9999.;
  fDeltaZ = -9999.;
  fRadius = -9999.;
  fOverlapY = -9999.;
  fOverlapZ = -9999.;
  fDeltaZ_Trigger = -9999.;

  //Begin temporary TODO: Delete later
  f2ndBestFlashT0 = -9999.;
  f2ndBestFlashTime = -9999.;
  f2ndBestFlashPEs = -9999.;
  f2ndBestFlashCenterY = -9999.;
  f2ndBestFlashCenterZ = -9999.;
  f2ndBestFlashWidthY = -9999.;
  f2ndBestFlashWidthZ = -9999.;
  f2ndBestDeltaT = -9999.;
  f2ndBestDeltaY = -9999.;
  f2ndBestDeltaZ = -9999.;
  f2ndBestRadius = -9999.;
  f2ndBestOverlapY = -9999.;
  f2ndBestOverlapZ = -9999.;
  //End temporary

} //End InitializeSlice()


double TPCPMTMatchingAna::CentroidOverlap(double center1, double center2, double width1, double width2) {
  //Centroid 2 is contained within Centroid 1, so overlap is the whole Centroid 2
  if ( (center1 - width1 < center2 - width2) && (center1 + width1 > center2 + width2) ) return (2 * width2);

  //Centroid 1 is contained within Centroid 2, so overlap is the whole Centroid 1
  else if ( (center1 - width1 > center2 - width2) && (center1 + width1 < center2 + width2) ) return (2 * width1);

  double difference = center1 - center2;
  if ( center1 > center2 ) difference *= -1;
  return difference + width1 + width2;
} //End CentroidOverlap()


double TPCPMTMatchingAna::CalculateAsymmetry(recob::OpFlash flash, int cryo) {
  double sumEast = 0.;
  double sumWest = 0.;

  //East Cryo flashes have a 180-element PE vector; 0-89 --> East wall, 90-179 --> West wall
  //West Cryo flashes have a 360-element PE vector; 0-179 --> ALL 0, 180-269 --> East wall, and 270-359 --> West wall
  int countingOffset = 0;
  if ( cryo == 1 ) countingOffset += 180;

  for ( int PMT = 0; PMT < 180; PMT++ ) {
    if ( PMT <= 89 ) sumEast += flash.PEs().at(PMT + countingOffset);
    else sumWest += flash.PEs().at(PMT + countingOffset);
  }

  return (sumWest - sumEast) / (sumWest + sumEast);
} //End CalculateAsymmetry


DEFINE_ART_MODULE(TPCPMTMatchingAna)
