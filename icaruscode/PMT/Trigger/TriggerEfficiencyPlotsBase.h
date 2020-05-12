/**
 * @file   icaruscode/PMT/Trigger/TriggerEfficiencyPlotsBase.h
 * @brief  Base class for _art_modules plotting trigger efficiencies.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 30, 2020
 */

#ifndef ICARUSCODE_PMT_TRIGGER_TRIGGEREFFICIENCYPLOTSBASE_H
#define ICARUSCODE_PMT_TRIGGER_TRIGGEREFFICIENCYPLOTSBASE_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // ADCCounts_t
#include "icaruscode/PMT/Trigger/Algorithms/details/TriggerInfo_t.h"
#include "icaruscode/PMT/Trigger/Utilities/PlotSandbox.h"
#include "icaruscode/PMT/Trigger/Data/MultiChannelOpticalTriggerGate.h"

// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // simulation_time, ...
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/Utilities/quantities/spacetime.h" // microseconds, ...
#include "lardataalg/Utilities/quantities/energy.h" // gigaelectronvolt, ...
#include "lardataalg/Utilities/intervals_fhicl.h" // microseconds from FHiCL
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/CoreUtils/get_elements.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::CryostatID

// framework libraries
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/ConsumesCollector.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"

// ROOT libraries
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TEfficiency.h"

// C/C++ standard libraries
#include <iosfwd> // std::ostream
#include <algorithm> // std::transform
#include <atomic>
#include <map>
#include <vector>
#include <array>
#include <string>
#include <optional>
#include <functional> // std::function<>, std::reference_wrapper<>
#include <memory> // std::unique_ptr
#include <utility> // std::forward(), std::pair<>, std::move()
#include <cstddef> // std::size_t


//------------------------------------------------------------------------------
namespace icarus::trigger {
  
  using namespace util::quantities::time_literals; // ""_ns ...
  
  class TriggerEfficiencyPlotsBase; 
  
} // icarus::trigger


//------------------------------------------------------------------------------
namespace icarus::trigger::details {
  
  struct EventInfo_t;
  
  std::ostream& operator<< (std::ostream& out, EventInfo_t const& info);
  
  struct TreeHolder;
  struct EventIDTree;
  struct PlotInfoTree;
  struct EventInfoTree;
    
} // namespace icarus::trigger::details


//------------------------------------------------------------------------------
/// Information about the event.
struct icarus::trigger::details::EventInfo_t {
  
  using GeV = util::quantities::gigaelectronvolt;
  
  /// Constructor. As if nobody noticed.
  EventInfo_t() { fInteractions.fill(0U); }
  
  // --- BEGIN Query interface -----------------------------------------------
  /// @name Query interface
  /// @{

  /// Returns the number of weak charged current interactions in the event.
  unsigned int nWeakChargedCurrentInteractions() const
    { return fInteractions[itWCC]; }

  /// Returns the number of weak neutral current interactions in the event.
  unsigned int nWeakNeutralCurrentInteractions() const
    { return fInteractions[itWNC]; }

  /// Returns the number of weak current interactions in the event.
  unsigned int nWeakCurrentInteractions() const
    {
      return
        nWeakChargedCurrentInteractions() + nWeakNeutralCurrentInteractions();
    }

  /// Returns whether the event is generated as a neutrino CC interaction.
  bool isWeakChargedCurrent() const
    { return nWeakChargedCurrentInteractions() > 0U; }

  /// Returns whether the event is generated as a neutrino NC interaction.
  bool isWeakNeutralCurrent() const
    { return nWeakNeutralCurrentInteractions() > 0U; }

  /// Returns whether the event is generated as a neutrino interaction.
  bool isNeutrino() const
    { return isWeakChargedCurrent() || isWeakNeutralCurrent(); }
  
  /// Returns which neutrino flavor is present in an event
  bool isNu_mu() const { return nu_mu; }
  bool isNu_e() const { return nu_e; }

  /// Returns the neutrino PDG code
  int NeutrinoPDG() const { return fNeutrinoPDG; }

  /// Returns the total energy deposited in the detector during the event [GeV]
  GeV DepositedEnergy() const { return fEnergyDepTotal; }
  
  /// Returns the total energy deposited in the detector during beam [GeV]
  GeV DepositedEnergyInSpill() const { return fEnergyDepSpill; }

  /// Returns the energy deposited in the active volume during the event [GeV]
  GeV DepositedEnergyInActiveVolume() const { return fEnergyDepActive; }
  
  /// Returns the energy deposited in the active volume during the beam [GeV]
  GeV DepositedEnergyInSpillInActiveVolume() const
    { return fEnergyDepSpillActive; }

  // Returns the neutrino energy [GeV]
  GeV NeutrinoEnergy() const { return fNeutrinoEnergy; }

  // Returns the lepton energy [GeV]
  GeV LeptonEnergy() const { return fLeptonEnergy; }

  // Returns the interaction type
  int InteractionType() const { return fInteractionType; }

  /// Returns whether this type of event has a known vertex.
  bool hasVertex() const { return !fVertices.empty(); }
  
  /// Returns the number of known interaction vertices.
  unsigned int nVertices() const { return fVertices.size(); }
  
  /// Returns the list of a known interaction vertex.
  std::vector<geo::Point_t> const& Vertices() const { return fVertices; }
  
  /// Returns whether there is an interaction within the active volume.
  bool isInActiveVolume() const { return fInActiveVolume; }
  
  /// @}
  // --- END Query interface -------------------------------------------------


  // --- BEGIN Set interface -------------------------------------------------
  /// @name Set interface
  /// @{

  /// Marks this event as including _n_ more weak charged current interactions.
  void AddWeakChargedCurrentInteractions(unsigned int n = 1U)
    { fInteractions[itWCC] += n; }

  /// Marks this event as including _n_ more weak neutral current interactions.
  void AddWeakNeutralCurrentInteractions(unsigned int n = 1U)
    { fInteractions[itWNC] += n; }

  /// Marks the neutrino flavor
  void SetNu_mu(bool numu) { nu_mu = numu; }
  void SetNu_e(bool nue) { nu_e = nue; }

  /// Marks this event's neutrino type
  void SetNeutrinoPDG(int NU) { fNeutrinoPDG = NU; }

  /// Sets the total deposited energy of the event [GeV]
  void SetDepositedEnergy(GeV e) { fEnergyDepTotal = e; }

  /// Sets the energy of the event deposited during beam gate [GeV]
  void SetDepositedEnergyInSpill(GeV e) { fEnergyDepSpill = e; }

  /// Sets the total deposited energy of the event in active volume [GeV]
  void SetDepositedEnergyInActiveVolume(GeV e) { fEnergyDepActive = e; }

  /// Sets the energy of the event deposited during beam gate in active volume
  /// [GeV]
  void SetDepositedEnergyInActiveVolumeInSpill(GeV e)
    { fEnergyDepSpillActive = e; }

  // Sets the neutrino energy.
  void SetNeutrinoEnergy(GeV eNu) { fNeutrinoEnergy = eNu; }

  // Sets the lepton energy.
  void SetLeptonEnergy(GeV eL) { fLeptonEnergy = eL; }

  // Sets the interaction type
  void SetInteractionType(int type) { fInteractionType = type; }

  /// Set whether the event has relevant activity in the active volume.
  void SetInActiveVolume(bool active = true) { fInActiveVolume = active; }
  
  /// Adds a point to the list of interaction vertices in the event.
  void AddVertex(geo::Point_t const& vertex) { fVertices.push_back(vertex); }
  
  /// @}
  // --- END Set interface ---------------------------------------------------

  /// Prints the content of the object into a stream.
  void dump(std::ostream& out) const;

    private:

  // --- BEGIN interaction type constants ------------------------------------
  
  static constexpr std::size_t itWCC { 0U }; ///< Charged weak current.
  static constexpr std::size_t itWNC { 1U }; ///< Neutral weak current.
  static constexpr std::size_t NInteractionTypes { 2U };

  // --- END interaction type constants --------------------------------------

  std::array<unsigned int, NInteractionTypes> fInteractions;

  GeV fEnergyDepTotal       { 0.0 }; ///< Total deposited energy.
  GeV fEnergyDepSpill       { 0.0 }; ///< Energy deposited in spill.
  GeV fEnergyDepActive      { 0.0 }; ///< Energy deposited in active volume.
  /// Energy deposited in active volume in spill.
  GeV fEnergyDepSpillActive { 0.0 };

  int fNeutrinoPDG { 0 };
  int fInteractionType { 0 };

  GeV fNeutrinoEnergy { 0.0 };
  GeV fLeptonEnergy { 0.0 };
  //GeV fNucleonEnergy { 0.0 };
  
  bool nu_mu { false };
  bool nu_e { false };
  
  /// Whether the event has activity inside the active volume.
  bool fInActiveVolume { false };
  
  std::vector<geo::Point_t> fVertices; ///< Position of all vertices.

}; // struct icarus::trigger::details::EventInfo_t

inline std::ostream& icarus::trigger::details::operator<<
  (std::ostream& out, EventInfo_t const& info)
  { info.dump(out); return out; }


// --- BEGIN -- ROOT tree helpers ----------------------------------------------
struct icarus::trigger::details::TreeHolder {

  TreeHolder() = default;
  TreeHolder(TTree& tree): fTree(&tree) {}

  TTree& tree() { return *fTree; }
  TTree const& tree() const { return *fTree; }

    private:
  TTree* fTree = nullptr;

}; // struct TreeHolder


/**
 * @brief Class managing the serialization of event ID in a simple ROOT tree.
 *
 * The tree is supplied by the caller.
 * This object will create the proper branches into the tree and assign
 * addresses to them. Then it will assume they will stay assigned.
 *
 * On `assignEvent()`, the branch addresses are assigned the values from the
 * event ID. The tree is not `Fill()`-ed.
 *
 * The tree structure is: `Run/i:SubRun/i:Event/i`, with a single branch per
 * element.
 */
struct icarus::trigger::details::EventIDTree: public TreeHolder {

  /// Creates the required branches and assigns addresses to them.
  EventIDTree(TTree& tree);

  /// Fills the information of the specified event.
  void assignID(art::EventID const& id);
  void assignEvent(art::Event const& event) { assignID(event.id()); }

  UInt_t fRunNo;
  UInt_t fSubRunNo;
  UInt_t fEventNo;

}; // struct icarus::trigger::details::EventIDTree


/**
 * @brief Class managing the serialization of plot information in a simple ROOT
 *        tree.
 *
 * The tree is supplied by the caller.
 * This object will create the proper branches into the tree and assign
 * addresses to them. Then it will assume they will stay assigned.
 *
 * On `assign()`, the branch addresses are assigned the values from the
 * arguments. The tree is not `Fill()`-ed.
 *
 * The tree structure is:
 * `InPlots/O`,
 * with a single branch per element.
 *
 * Branches:
 *  * `InPlots` (bool): the event was _not_ filtered out before plotting
 *    (it may still belong to no category and eventually appear in no plot)
 *
 */
struct icarus::trigger::details::PlotInfoTree: public TreeHolder {

  /// Creates the required branches and assigns addresses to them.
  PlotInfoTree(TTree& tree);

  /**
   * @brief Fills the information of the specified event.
   * @param inPlots whether this event is plotted (as opposed to filtered out)
   */
  void assign(bool inPlots);

  Bool_t fInPlots;

}; // struct icarus::trigger::details::PlotInfoTree


/**
 * @brief Class managing the serialization of event information in a simple ROOT
 *        tree.
 *
 * The tree is supplied by the caller.
 * This object will create the proper branches into the tree and assign
 * addresses to them. Then it will assume they will stay assigned.
 *
 * On `assignEvent()`, the branch addresses are assigned the values from the
 * event information. The tree is not `Fill()`-ed.
 *
 * The tree structure is:
 * `CC/i:NC/i:IntType:I/NuE:D/OutLeptE:D/TotE/D:SpillE/D:InActive/O`,
 * with a single branch per element.
 *
 * Branches:
 *  * `CC` (unsigned integer): number of neutrino CC interactions in the event
 *  * `NC` (unsigned integer): number of neutrino NC interactions in the event
 *  * `IntType` (integer): code of interaction type (see `simb::int_type_`)
 *  * `NuE` (double): energy of the generated initial state neutrino [GeV]
 *  * `OutLeptE` (double): energy of the generated final state lepton [GeV]
 *  * `TotE` (double): total deposited energy in the event [GeV]
 *  * `SpillE` (double): total deposited energy during the beam gate [GeV]
 *  * `InActive` (bool): whether an interaction happened in active volume'
 *      this requires an interaction vertex (e.g. cosmic rays are out)
 *  * `NVertices` (unsigned integer): number of interaction vertices in event
 *  * `Vertices_` (list of points): the location of all the interaction vertices
 *    in the event; it's a vector of GenVector 3D points (can access coordinates
 *    as `Vertices.X()` or `Vertices.fCoordinates.fX`)
 *
 */
struct icarus::trigger::details::EventInfoTree: public TreeHolder {

  /// Creates the required branches and assigns addresses to them.
  EventInfoTree(TTree& tree);

  /**
   * @brief Fills the information of the specified event.
   * @param info event information to fill the tree with
   * @param inPlots whether this event is plotted (as opposed to filtered out)
   */
  void assignEvent(EventInfo_t const& info);

  UInt_t fCC;
  UInt_t fNC;
  Int_t fIntType;
  Double_t fNuE;
  Double_t fOutLeptE;
  Double_t fTotE;
  Double_t fSpillE;
  Double_t fActiveE;
  Double_t fSpillActiveE;
  UInt_t fNVertices;
  std::vector<geo::Point_t> fVertices; // is this ROOT tree friendly?
  
  Bool_t fInActive;
  
}; // struct icarus::trigger::details::EventInfoTree


// --- END -- ROOT tree helpers ------------------------------------------------


//------------------------------------------------------------------------------
/**
 * @brief Helper class to produce plots about trigger simulation and trigger
 *        efficiency.
 * @see icarus::trigger::MajorityTriggerEfficiencyPlots
 * 
 * This helper class provides the foundation for writing a module producing sets
 * of plots based on trigger primitives given in input.
 * 
 * The following documentation mainly deals with the most standard configuration
 * and operations in ICARUS. Modules made with this helper and the modules
 * upstream of them in data processing have quite some knobs that can be
 * manipulated to test unorthodox configurations.
 * 
 * An example of module implemented with this class is
 * `icarus::trigger::MajorityTriggerEfficiencyPlots`.
 * 
 * 
 * Overview of trigger simulation steps
 * =====================================
 * 
 * This simulation only trigger primitives derived from the optical detector
 * via hardware (V1730B boards).
 * The trigger simulation branches from the standard simulation of optical
 * detector waveforms (`icarus::SimPMTIcarus` module).
 * From there, multiple steps are needed.
 * 
 * 1. Produce single-PMT-channel discriminated waveforms: the discrimination
 *    threshold can be chosen ("ADC threshold" or just "threshold" in the
 *    following), and the result is one binary discriminated waveform that
 *    has value `0` when the waveform is under threshold and `1` when it's
 *    above threshold, with the same time discretization as the PMT waveform.
 *    Each discriminated waveform spans the whole event readout time and
 *    therefore it may merge information from multiple readout waveforms
 *    (`raw::OpDetWaveform`), but all from the same optical detector channel.
 *    This step can be performed with the module
 *    `icarus::trigger::DiscriminatePMTwaveforms`.
 * 2. Combine the discriminated waveforms in pairs, in the same fashion as the
 *    V1730B readout board does to produce LVDS output. The output of this step
 *    is roughly half as many discriminated outputs as there were from the
 *    previous step (some channels are not paired). This output will be called
 *    _trigger primitives_ because it is what the trigger logics is based on.
 *    We say that a trigger primitive is "on" when its level is `1`.
 *    This step can be performed with the module `icarus::trigger::LVDSgates`.
 * 3. Simulate the trigger logic based on the trigger primitives _(see below)_.
 *    This usually includes the requirement of coincidence with the beam gate.
 * 
 * Trigger logic may be complex, being implemented in a FPGA.
 * Many options are available, including:
 * 
 * * coincidence of a minimum number of trigger primitives: the event is trigger
 *   if at least _N_ trigger primitives are on at the same time;
 * * sliding window: primitives are grouped depending on their location, and
 *   there is a requirement on the minimum number of primitives in "on" level
 *   within one or more of these groups. For example, any sliding window of 30
 *   channels (corresponding to 16 trigger primitives in standard ICARUS PMT
 *   readout installation) has to contain at least 10 trigger primitives "on"
 *   at the same time; or there must be two sliding windows with that condition;
 *   etc. Sliding windows can be obtained from further processing of the LVDS
 *   trigger primitives, for example with the module
 *   `icarus::trigger::SlidingWindowTrigger`.
 * 
 * This class is expected to support modules implementing different trigger
 * logics.
 * 
 * 
 * Data objects for discriminated waveforms
 * -----------------------------------------
 * 
 * @anchor TriggerEfficiencyPlotsBase_Data
 *
 * A discriminated waveform is the information whether the level of a waveform
 * is beyond threshold, as function of time.
 * A discriminated waveform may be binary, i.e. with only levels `0` and `1`
 * based on a single threshold, or with multiple levels.
 * Also the numerical _addition_ of two binary discriminated waveforms
 * yields a multi-level waveform (in fact, three levels -- `0`, `1` and `2`).
 * 
 * We represent this data in the forms of "events": an event is a change of
 * level happening at a certain time. The class holding this information,
 * `icarus::trigger::TriggerGateData`, covers the whole time, starting with a
 * ground level `0`. The next event will be a "opening" that increases the
 * level, usually to `1`. Other changing events may follow, and typically the
 * last one will bring the level back to `0`.
 * 
 * This information is joined by a list of _channel numbers_ in order to
 * represent a discriminated waveform e.g. from the optical detector.
 * There may be one or more channels associated to a discriminated waveform,
 * but for us there should always be at least one.
 * The discriminated waveforms from PMT readout (`raw::OpDetWaveform`) are
 * associated to a single channel, while LVDS trigger primitives are most often
 * associated to two channels (some channels are not paired and will have only
 * one channel). A sliding window will have as many channels as the PMT it
 * covers. The global trigger should have _all_ channels, while in ICARUS each
 * of the two discriminated wabeforms from a cryostat should have half the
 * channels.
 * This information is represented in the class
 * `icarus::trigger::ReadoutTriggerGate`, which inherits from
 * `icarus::trigger::TriggerGateData`.
 * This class is generic and can hold any representation for the time of the
 * level changing events, for the levels, and for the identifiers of the
 * channels. ICARUS specifies a data type for each of these quantities, and
 * the resulting `icarus::trigger::ReadoutTriggerGate` class instance is called
 * `icarus::trigger::OpticalTriggerGateData_t`.
 * 
 * @note The class `icarus::trigger::OpticalTriggerGateData_t` is the one that
 *       gets written to disk in the _art_ ROOT files.
 *       That is not a class by itself, but rather an alias of
 *       `icarus::trigger::ReadoutTriggerGate`, and many coding tools will call
 *       it in the latter way.
 * 
 * The class `icarus::trigger::OpticalTriggerGate` is currently the most
 * commonly used in the code. It adds to the information of
 * `icarus::trigger::ReadoutTriggerGate`, from which it derives, a list of
 * optical waveforms (`raw::OpDetWaveform`) it originates from.
 * 
 * Finally, the classes `icarus::trigger::SingleChannelOpticalTriggerGate` and
 * `icarus::trigger::MultiChannelOpticalTriggerGate` do not add any information
 * to the `icarus::trigger::OpticalTriggerGate` they derive from, but they
 * have an interface explicitly tuned for discriminated waveforms from a single
 * channel or from multiple channels, respectively (for example, the former
 * provides a `channel()` method returning a single channel identifier, while
 * the latter provides a `channels()` method returning a list of channels).
 * 
 * These three higher level classes, `icarus::trigger::OpticalTriggerGate` and
 * derivatives, _can't be directly saved_ in _art_ ROOT files.
 * There are utilities available in
 * `icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h` that can convert them
 * into a collection of `icarus::trigger::OpticalTriggerGateData_t` objects
 * plus a collection of _art_ associations (for writing), and the other way
 * around (for reading). The module `icarus::trigger::LVDSgates` uses both sides
 * and can be used as an illustration of the functionality.
 * 
 * A module is provided, called `icarus::trigger::DumpTriggerGateData`, which
 * dumps on screen or on text file the information from a collection of
 * discriminated waveforms in a (sort-of) human readable format.
 * An example configuration for this module is provided in `icaruscode`, called
 * `dump_triggergatedata_icarus.fcl`.
 * 
 * The main functions to manipulate the trigger gates are defined in the very
 * base class, `icarus::trigger::TriggerGateData`: these allow e.g. to find
 * events and query the level at a given time.
 * Another important set of features is also present in
 * `icarus::trigger::TriggerGateData` and replicated in the higher levels:
 * the combination of trigger gates by sum (`OR`), multiplication (_AND_),
 * minimum and maximum value.
 * 
 * @note Combining a multi-level gate via `Min()` with a binary gate results
 *       into a binary gate which is logic AND of the two gates.
 * 
 * 
 * On the terminology
 * -------------------
 * 
 * In this documentation, the "discriminated waveforms" are sometimes called
 * "trigger gates" and sometimes "trigger primitives".
 * Although these concepts are, strictly, different, usually the difference
 * does not affect the meaning and they end up being exchanged carelessly.
 * 
 * 
 * Trigger logic algorithm
 * ========================
 * 
 * @anchor TriggerEfficiencyPlotsBase_Algorithm
 * 
 * The trigger logic algorithm is completely absent from this helper, and it is
 * entirely delegated to the derived classes (`simulateAndPlot()`). There are
 * nevertheless a lot of functions that may help to process the input trigger
 * primitives and fill some standard plots. Many of them are listed in the
 * documentation of `simulateAndPlot()`.
 * 
 * 
 * Event categories
 * =================
 * 
 * Each event is assigned to a set of categories, and its information
 * contributes to the plots of all those categories.
 * The categories are passed to the helper in the `initializePlots()` call to
 * create all the plots to be filled later.
 * A default list of categories is provided as `DefaultPlotCategories`.
 * 
 * 
 * Module usage
 * =============
 * 
 * The following aspects are common to all modules based on this helper.
 * 
 * 
 * Input data products
 * --------------------
 * 
 * This module uses the following data products as input:
 * 
 * * trigger primitives:
 *     * `std::vector<icarus::trigger::OpticalTriggerGateData_t>` (labels out of
 *       `TriggerGatesTag` and `Thresholds`): full sets of discriminated
 *       waveforms, each waveform possibly covering multiple optical channels,
 *       and their associations to optical waveforms. One set per threshold;
 *       note that these are converted into `InputTriggerGate_t` data type for
 *       internal use;
 *     * associations with `raw::OpDetWaveform` (currently superfluous);
 * * event characterization:
 *     * `std::vector<simb::MCTruth>`: generator information (from
 *       `GeneratorTags`; multiple generators are allowed at the same time);
 *       it is used mostly to categorize the type of event (background,
 *       weak charged current, electron neutrino, etc.);
 *     * `std::vector<simb::MCParticle>`: particles propagating in the detector
 *       (from `PropagatedParticles`); currently not used;
 *     * `std::vector<sim::SimEnergyDeposit>`: energy deposited in the active
 *       liquid argon volume (from `EnergyDepositTags`); it is used to
 *       quantify the energy available to be detected in the event.
 * 
 * 
 * Output plots
 * -------------
 * 
 * @anchor TriggerEfficiencyPlotsBase_Plots
 * 
 * The module produces a standard set of plots for each configured ADC threshold
 * and for each event category.
 * The plots are saved via _art_ service `TFileService` in a ROOT file and they
 * are organized in nested ROOT directories under the module label directory
 * which is assigned by _art_:
 * 
 * * `Thr###` (outer level) describes the ADC threshold on the discriminated
 *   waveforms: the threshold is `###` from the baseline;
 * * `\<Category\>` (inner level) describes the category of events included in
 *   the plots (e.g. `All`, `NuCC`; see `PlotCategories`).
 * 
 * Each of the inner ROOT directories contains a full set of plots, whose name
 * is the standard plot name followed by its event category and threshold
 * (e.g. `Eff_NuCC_Thr15` for the trigger efficiency plot on neutrino charged
 * current events with 15 ADC counts as PMT threshold).
 * 
 * Each set of plots, defined in
 * `icarus::trigger::TriggerEfficiencyPlots::initializePlotSet()`
 * (and, within, `initializeEventPlots()`
 * and `initializeEfficiencyPerTriggerPlots()`), is contributed only by the
 * events in the set category.
 * 
 * There are different "types" of plots. Some
 * @ref TriggerEfficiencyPlotsBase_SelectionPlots "do not depend on triggering at all",
 * like the deposited energy distribution. Others
 * @ref TriggerEfficiencyPlotsBase_MultiTriggerPlots "cross different trigger definitions",
 * like the trigger efficiency as function of trigger requirement. Others still
 * @ref TriggerEfficiencyPlotsBase_SingleTriggerPlots "assume a single trigger definition":
 * this is the case of trigger efficiency plots versus energy. Finally, there are
 * @ref TriggerEfficiencyPlotsBase_SingleTriggerResponsePlots "plots that depend on a specific trigger definition and outcome":
 * this is the case of all the plots including only triggering or non-triggering
 * events.
 * 
 * A list of plots follows for each plot type.
 * All the plots are always relative to a specific optical detector channel
 * threshold (ADC) and a broad event category.
 * 
 * 
 * ### Plots independent of the triggers (selection plots)
 * 
 * @anchor TriggerEfficiencyPlotsBase_SelectionPlots
 * 
 * These plots are stored directly in a threshold/category folder:
 * 
 * * `EnergyInSpill`: total energy deposited in the detector during the time the
 *   beam gate is open. It is proportional to the amount of scintillation light
 *   in the event;
 * * `EnergyInSpillActive`: energy deposited in the active volume of the
 *   detector during the time the beam gate is open; the active volume is
 *   defined as the union of all TPC drift voulmes;
 * * plots specific to neutrino interactions (if the event is not a neutrino
 *   interaction, it will not contribute to them); if not otherwise specified,
 *   only the first neutrino interaction in the event is considered:
 *     * `InteractionType`: code of the interaction type, as in
 *       `sim::TruthInteractionTypeName`;
 *     * `NeutrinoEnergy`: generated energy of the interacting neutrino;
 *     * `LeptonEnergy`: generated energy of the lepton out of the _first_
 *       neutrino interaction;
 *     * `InteractionVertexYZ`: coordinates of the location of all interactions
 *       in the event, in world coordinates, projected on the anode plane.
 * 
 * @note In fact, these plots usually do not even depend on the ADC threshold
 *       of the optical detector channels. Nevertheless, they are stored in the
 *       folders under specific thresholds, and they are duplicate.
 * 
 * 
 * ### Plots including different trigger requirements
 * 
 * @anchor TriggerEfficiencyPlotsBase_MultiTriggerPlots
 * 
 * These plots collect information from scenarios with different trigger
 * requirements, but still with the same underlying optical detector channel
 * ADC threshold.
 * No plots are provided in this category by `TriggerEfficiencyPlotsBase` at
 * the moment. Implementing modules can add plots by overriding the method
 * `initializePlotSet()`:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * void MyTriggerEfficiencyPlots::initializePlotSet
 *   (PlotSandbox& plots, std::vector<SettingsInfo_t> const& settings) const
 * {
 *   // (inherited)
 *   helper().TriggerEfficiencyPlotsBase::initializePlotSet(plots, settings);
 *   
 *   // ... add more plot definitions, e.g.
 *   plots.make<TEfficiency>(
 *     "Eff",
 *     "Efficiency of triggering;requested primitives;trigger efficiency",
 *     settings.size(), 0, settings.size()
 *     );
 *   
 * } // MyTriggerEfficiencyPlots::initializePlotSet()
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * 
 * ### Plots depending on a specific trigger definition
 * 
 * @anchor TriggerEfficiencyPlotsBase_SingleTriggerPlots
 * 
 * The following plots depend on the full definition of the trigger, including
 * PMT thresholds _and_ other requirements.
 * They are hosted each in a subfolder of the threshold/category folder, with
 * a name encoding the requirement: `ReqXX` for trigger with minimum required
 * primitives `XX`.
 * 
 * * `EffVsEnergyInSpill`: trigger efficiency as function of the total energy
 *   deposited during the beam gate;
 * * `EffVsEnergyInSpillActive`: trigger efficiency as function of the energy
 *   deposited in the TPC active volumes during the beam gate;
 * * `EffVsNeutrinoEnergy`, `EffVsLeptonEnergy`: trigger efficiency as function
 *   of the true energy of the first interacting neutrino, and of the outgoing
 *   lepton in the final state of the interaction, respectively;
 * * `TriggerTick`: time of the earliest trigger. Only triggering events
 *   contribute.
 * 
 * The parameters are defined in the same way as in the
 * @ref TriggerEfficiencyPlotsBase_SelectionPlots "selection plots", unless stated
 * otherwise.
 *
 *
 * ### Plots depending on a specific trigger definition and response
 *
 * @anchor TriggerEfficiencyPlotsBase_SingleTriggerResponsePlots
 *
 * These plots depend on the trigger definition, as the ones in the previous
 * type, and on the actual trigger response.
 * They are hosted each in a subfolder of the threshold/category/requirement
 * folder, with a name encoding the response: `triggering` for triggering
 * events, `nontriggering` for the others.
 * 
 * Their event pool is filtered to include only the events in the current
 * category which also pass, or do not pass, the trigger requirement.
 * 
 * The plots are:
 * * `EnergyInSpill`, `EnergyInSpillActive`, `InteractionType`,
 *   `NeutrinoEnergy`,`LeptonEnergy`, `InteractionVertexYZ`:
 *   these are defined in the same way as the
 *   @ref TriggerEfficiencyPlotsBase_SelectionPlots "selection plots"
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * @anchor TriggerEfficiencyPlotsBase_Configuration
 * 
 * A terse description of the parameters is printed by running
 * `lar --print-description TriggerEfficiencyPlots`.
 * 
 * * `TriggerGatesTag` (string, mandatory): name of the module
 *     instance which produced the trigger primitives to be used as input;
 *     it must not include any instance name, as the instance names will be
 *     automatically added from `Thresholds` parameter.
 *     The typical trigger primitives used as input may be LVDS discriminated
 *     output (e.g. from `icarus::trigger::LVDSgates` module) or combinations
 *     of them (e.g. from `icarus::trigger::SlidingWindowTrigger` module).
 * * `Thresholds` (list of integers, mandatory): list of the discrimination
 *     thresholds to consider, in ADC counts. A data product containing a
 *     digital signal is read for each one of the thresholds, and the tag of the
 *     data product is expected to be the module label `TriggerGatesTag` with as
 *     instance name the value of the threshold (e.g. for a threshold of 60 ADC
 *     counts the data product tag might be `LVDSgates:60`).
 * * `GeneratorTags` (list of input tags, default: `[ generator ]`): a list of
 *     data products containing the particles generated by event generators;
 * * `DetectorParticleTag` (input tag, default: `largeant`): data product
 *     containing the list of particles going through the argon in the detector;
 * * `EnergyDepositTags`
 *     (list of input tags, default: `[ "largeant:TPCActive" ]`): a list of
 *     data products with energy depositions;
 * * `BeamGateDuration` (time, _mandatory_): the duration of the beam
 *     gate; _the time requires the unit to be explicitly specified_: use
 *     `"1.6 us"` for BNB, `9.5 us` for NuMI (also available as
 *     `BNB_settings.spill_duration` and `NuMI_settings.spill_duration` in
 *     `trigger_icarus.fcl`);
 * * `TriggerTimeResolution` (time, default: `8 ns`): time resolution for the
 *     trigger primitives;
 * * `EventTreeName` (optional string): if specified, a simple ROOT tree is
 *     created with the information from each event (see `EventInfoTree` for
 *     its structure);
 * * `EventDetailsLogCategory` (optional string): if specified, information for
 *     each single event is output into the specified stream category; if
 *     the string is specified empty, the default module stream is used, as
 *     determined by `LogCategory` parameter;
 * * `LogCategory` (string, default `TriggerEfficiencyPlots`): name of category
 *     used to stream messages from this module into message facility.
 * 
 * An example job configuration is provided as `maketriggerplots_icarus.fcl`.
 * 
 * 
 * Technical description of the module
 * ====================================
 * 
 * The modules based on this helper read
 * @ref TriggerEfficiencyPlotsBase_Data "trigger gate data products" from
 * different ADC thresholds, and for each threahold they combine the data into a
 * trigger response depending on a ADC threshold.
 * Then they apply different trigger settings and for each one they fill plots.
 *
 * All the input sets (each with its own ADC treshold) are treated independently
 * from each other.
 * 
 * Each event is assigned to event categories. These categories do not depend on
 * trigger definition nor primitives. For example, a electron neutrino neutral
 * current interaction event would be assigned to the neutral current neutrino
 * interaction category, to the electron neutrino category, and so on.
 * A set of plot is produced at once for each of the event categories.
 *
 * This module is designed as follows:
 *
 * * each _art_ event is treated independently from the others (this is the
 *   norm; method: `process()`, to be called in module's `analyze()`);
 * * information for the plots which do not depend on the trigger primitives
 *   are extracted once per event (method: `extractEventInfo()`);
 * * the event is assigned its categories (method: `selectedPlotCategories()`,
 *   with the categories defined at the time of `initializePlots()` call);
 * * for each input primitive set (i.e. for each ADC threshold):
 *     * trigger logic is applied, the beam gate is imposed, a trigger result
 *       is extracted and plots are filled; all is done in the method
 *       `simulateAndPlot()` provided by the derived module.
 *
 * Note that the plots depending on the response may require multiple fillings.
 * For example, a trigger efficiency plot as function of the trigger settings
 * is a single plot in the plot set which requires a filling for every
 * trigger setting.
 * In a similar fashion, the plot of trigger time requires one filling for
 * each requirement level that is met.
 *
 *
 * Organization of the plots
 * --------------------------
 * 
 * @anchor TriggerEfficiencyPlotsBase_OrganizationOfPlots
 *
 * Plots are written on disk via the standard _art_ service `TFileService`,
 * which puts them in a ROOT subdirectory named as the instance of this module
 * being run.
 *
 * As described above, there are two "dimensions" for the plot sets: there is
 * one plot set for each ADC threshold and for each event category, the total
 * number of plot sets being the product of the options in the two dimensions.
 *
 * The code is structured to work by threshold by threshold, because the trigger
 * response depends on the threshold but not by the event category: the response
 * is computed once per threshold, then all the plots related to that response
 * (including all the event categories) are filled.
 *
 * The structure in the file reflects this logic, and there are two levels of
 * ROOT directories inside the one assigned by _art_ to this module:
 * * the outer level pins down the ADC threshold of the plots inside it;
 *   the name of the directory follows the pattern `Thr###` (### being the ADC
 *   threshold), which is the "threshold tag";
 * * the inner level is the category of events included in the plots, and the
 *   name of the directories reflect the name of the category as defined in
 *   the corresponding `PlotCategory` object (`PlotCategory::name()`); this
 *   defines the "category tag".
 *
 * In each inner directory, a complete set of plots is contained.
 * The name of each plot is a base name for that plot (e.g. `Eff` for the
 * efficiency plot) with the event category tag and the threshold tag appended
 * (separated by an underscore character, `"_"`). The title of the plot is also
 * modified to include the description of the plot category (as defined in
 * `PlotCategory::description()`) and the ADC threshold.
 * Therefore, within the same module directory, all plot names are different.
 *
 *
 * Adding a plot
 * --------------
 *
 * When adding a plot, two actions are necessary:
 *
 * 1. initialize the new plot
 * 2. fill the new plot
 *
 * The initialization happens in
 * `icarus::trigger::TriggerEfficiencyPlots::initializePlotSet()` method.
 * A request must be issued to the
 * @ref TriggerEfficiencyPlotsBase_PlotSandboxes "plot sandbox" to "make" the
 * plot.
 * In general it can be expected that all the arguments `args` in a call
 * `plots.make<PlotType>(args)` are forwarded to the constructor of `PlotType`,
 * so that to make a new `PlotType` object one has to consult only the
 * constructor of that class. Be aware though that the library expects the first
 * argument to be the ROOT name of the new object and the second to be its ROOT
 * title.
 * It is possible to work around this requirement with additional coding.
 *
 * The method performing the plotting is `simulateAndPlot()`.
 * The plots should be filled inside loops going through all plot sets.
 * The existing code already does that in two loops, one for the plots
 * depending on the trigger response requirement, and the other for the plots
 * _not_ depending on it.
 *
 *
 * ### About plot sandboxes
 * 
 * @anchor TriggerEfficiencyPlotsBase_PlotSandboxes
 *
 * For the sake of this module, a plot sandbox is an object similar to a ROOT
 * `TDirectory`, which can mangle the objects it contains to change ("process")
 * their name and title, and that interacts properly with `art::TFileDirectory`
 * so make `TFileService` happy.
 * The processing of the name and title is used to add the category and
 * threshold tags to the plot ROOT name and the corresponding annotations to
 * the plot title, as described
 * @ref TriggerEfficiencyPlotsBase_OrganizationOfPlots "above".
 *
 *
 * Adding an event category
 * -------------------------
 * 
 * @anchor TriggerEfficiencyPlotsBase_AddingCategory
 *
 * Event categories are specified to `initializePlots()` call: each category is
 * described by a simple object `PlotCategory` which contains a name (used in
 * ROOT directories, plot name tags and to univocally identify the category),
 * a description (used in plot titles) and a condition for the event to fulfill
 * in order to qualify for the category. The condition is a function object
 * that accepts an `EventInfo_t` object with the relevant event information,
 * and returns its decision (`true` if the event belongs to this category).
 *
 * To add a category, it is enough to add an entry in the collection of plot
 * categories passed to `initializePlots()`. The `DefaultPlotCategories` list
 * can be used as a starting point.
 *
 * The condition function adjudicates based on the information in `EventInfo_t`.
 * It is likely though that if a new category is needed, the information
 * currently available in `EventInfo_t` is not enough to decide the response.
 * In that case, `EventInfo_t` must be extended to hold the newly required
 * information. In addition, the method `extractEventInfo()` must be modified
 * to set the new information. The documentation of that method expands on this
 * topic.
 *
 *
 * Usage of the helper
 * ====================
 * 
 * An example of how to turn this helper class in a full blown _art_ module is
 * `icarus::trigger::MajorityTriggerEfficiencyPlots`.
 *
 *
 * Changing the trigger logic
 * ---------------------------
 *
 * The trigger logic is required to be coded in `simulateAndPlot()` of the
 * derived class. Some details on implementing it can be found in the
 * documentation of the method itself.
 *
 * To implement a trigger logic, several components need to be coordinated:
 *
 * 1. combination of trigger primitives;
 * 2. the data structure the combination is stored in;
 * 3. interpretation of that data structure in terms of trigger firing
 *    (which we call here "trigger settings");
 * 4. configuration of the trigger logic from user parameters.
 *
 * A module may simulate a trigger which counts the number of primitives that
 * are "on" and imposing a required minimum on that number for the trigger to
 * fire. This is a possible implementation including those elements:
 *
 * 1. combination of trigger primitives: primitives are summed into a single
 *    multilevel gate (with beam gate on top);
 * 2. the data structure the combination is stored in: a
 *    @ref TriggerEfficiencyPlotsBase_Data "trigger gate object";
 * 3. interpretation of that data structure in terms of trigger firing:
 *    the logic requires a minimum number of primitives being on at the same
 *    time; the interpretation uses the level of the combined trigger gate as
 *    the number of primitives on at every time, and decides that the first
 *    time this number meets the requirement, the trigger fires;
 * 4. configuration of the trigger logic from user parameters:
 *    it's a parameter of the module itself (`MinimumPrimitives`).
 *
 * An example of how a different trigger logic could be implemented following
 * a similar model: let's assume we want a trigger based on two sliding windows.
 * A window is a group of LVDS signals from neighboring channels. We can split
 * the whole optical detector channels into many windows, and we can make the
 * windows overlap: for example, a window includes the first 30 channels behind
 * one anode plane, the second window includes 30 channels starting from the
 * fifteenth, the third window includes 30 channels starting from the thirtieth,
 * and so on. In ICARUS, 30 channels are served in 16 LVDS discriminated
 * waveforms, i.e. 16 trigger primitives. We may require as a loose trigger
 * that any window contains at least 10 primitives on at the same time (that is
 * about 20 PMT channels beyond ADC threshold), and as a tighter trigger that
 * there is a window with at least 8 trigger primitives on, _and_ the
 * corresponding window in the opposite TPC has at least 6 primitives on
 * at the same time. Note that a splitting 7/7 won't do: one of the two windows
 * must have 8 or more primitives on.
 * The input primitives from the point of view of the module now become the
 * "sliding window" combination of the LVDS trigger primitives, as produced e.g.
 * by `icarus::trigger::SlidingWindowTrigger` module.
 * How this can be implemented following the pattern of this module:
 *
 * 1. combination of trigger primitives: we prepare two primitives:
 *    one describes the maximum level among all the windows, and the other
 *    describes the level of the window opposite to the maximum (the way to
 *    achieve this is quite complicate);
 * 2. the data structure the combination is stored in: two
 *    @ref TriggerEfficiencyPlotsBase_Data "trigger gate objects";
 * 3. interpretation of that data structure in terms of trigger firing:
 *    for the loose trigger, we look for the first time the maximum window level
 *    reaches the designated threshold; for the tighter trigger, some
 *    additional logic is required, discriminating the maximum window level
 *    primitive to the higher threshold, the opposite window level on the lower
 *    threshold, and setting them in coincidence; then the time the resulting
 *    gate opens, if any, is the trigger time;
 * 4. configuration of trigger logic parameters: the two requirements (single
 *    window, double window with their minimum number of LVDS primitives)
 *    are read from the module FHiCL configuration; for computing efficiency,
 *    either in the constructor or in a `beginJob()` method we may also prepare
 *    the associations (pairing) between opposite windows or just channels;
 *
 * @todo add ability to discriminate a trigger gate object?
 *       discrimination on level _N_ (>=N -> 1, \<N -> 0) can be obtained with
 *       adding a _-N_ uniform level, flooring on level 0 and (if needed) maxing
 *       on level 1.
 *
 * 
 * Code style: quantities with units
 * ==================================
 * 
 * To avoid issues with missing or wrong conversions, this code often uses
 * LArSoft quantities library. A variable with a `Quantity` type is represented
 * with a specific unit (e.g. microseconds) and can undergo only limited
 * manipulation. The allowed manipulation should guarantee that the unit of
 * the quantity is correctly preserved. For example, it is not possible to
 * add to a `microseconds` interval a pure number (e.g. `9.0`), but rather
 * it is only possible to add time quantities (e.g. another `microseconds`
 * variable, or a `nanoseconds` variable, or a literal value with unit, like
 * `9.0_ns`), and those quantities are properly converted, with the caveat that
 * rounding errors may happen that a more careful explicit handling could avoid.
 * Also, there are two types of variables that can feature the same unit,
 * intervals and points. A point can't be scaled (e.g. you can't "double" the
 * time the beam gate opens, while you can double the beam gate _duration_)
 * and it can't be added to another point (the difference between two points
 * is an interval).
 * 
 * To avoid mistakes in the conversion between different time scales, the
 * LArSoft utility library `detinfo::DetectorTimings` is used, which is a
 * wrapper of `DetectorClocks` service provider that makes the interface to
 * convert times easier and uniform. Here we use it to convert time points
 * from the simulation (which are expressed in nanoseconds and are starting
 * at trigger time) into optical detector readout ticks, and vice versa.
 * The values returned by this library have encoded in them which time scale
 * they belong to, and in which unit they are measured.
 * 
 */
class icarus::trigger::TriggerEfficiencyPlotsBase {

  // no serviceable part for anyone except derived classes
  
    protected:
  
  using microseconds = util::quantities::intervals::microseconds;
  using nanoseconds = util::quantities::intervals::nanoseconds;

  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Sequence<art::InputTag> GeneratorTags {
      Name("GeneratorTags"),
      Comment("labels of the event generators"),
      std::vector<art::InputTag>{ "generator" }
      };

    fhicl::Atom<art::InputTag> DetectorParticleTag {
      Name("DetectorParticleTag"),
      Comment("label of simulated particles through the detector"),
      "largeant" // tradition demands
      };

    fhicl::Sequence<art::InputTag> EnergyDepositTags {
      Name("EnergyDeposits"),
      Comment("label of energy deposition data product(s) in the detector"),
      std::vector<art::InputTag>{ "largeant:TPCActive" }
      };

    fhicl::Atom<std::string> TriggerGatesTag {
      Name("TriggerGatesTag"),
      Comment("label of the input trigger gate data product (no instance name)")
      };

    fhicl::Sequence<raw::ADC_Count_t> Thresholds {
      Name("Thresholds"),
      Comment("thresholds to consider [ADC counts]")
      };

    fhicl::Atom<microseconds> BeamGateDuration {
      Name("BeamGateDuration"),
      Comment("length of time interval when optical triggers are accepted")
      };

    fhicl::Atom<nanoseconds> TriggerTimeResolution {
      Name("TriggerTimeResolution"),
      Comment("resolution of trigger in time"),
      8_ns
      };
    
    fhicl::Atom<bool> PlotOnlyActiveVolume {
      Name("PlotOnlyActiveVolume"),
      Comment
        ("only events within TPC active volume are plot (if that makes sense)"),
      true
      };
    
    fhicl::OptionalAtom<std::string> EventTreeName {
      Name("EventTreeName"),
      Comment("name of a ROOT tree where to store event-by-event information")
      };

    fhicl::OptionalAtom<std::string> EventDetailsLogCategory {
      Name("EventDetailsLogCategory"),
      Comment("name of the category used for event information output")
      };

    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment("name of the category used for the output"),
      "SlidingWindowTrigger" // default
      };

  }; // struct Config

  // --- END Configuration -----------------------------------------------------


  using EventInfo_t = details::EventInfo_t; // type alias
  using TriggerInfo_t = details::TriggerInfo_t; // type alias
  
  
  // --- BEGIN Constructors ----------------------------------------------------

  /// Constructor; requires a configuration and module's `consumesCollector()`.
  explicit TriggerEfficiencyPlotsBase
    (Config const& config, art::ConsumesCollector& consumer);

  // --- END Constructors ------------------------------------------------------


  // --- BEGIN Framework hooks -------------------------------------------------

  /// Fills the plots. Also extracts the information to fill them with.
  void process(art::Event const& event);
  
  /// Prints end-of-job summaries.
  void printSummary() const;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
  // --- BEGIN Helper interface ------------------------------------------------
  //
  // Mediated access to private stuff.
  //
  
  /// Returns the name of the log category.
  std::string const& logCategory() const { return fLogCategory; }
  
  /// Returns a pointer to the tree where event information is written.
  TTree* eventTree() { return fIDTree? &(fIDTree->tree()): nullptr; }
  
  /// Returns the number of configured PMT thresholds.
  std::size_t nADCthresholds() const { return fADCthresholds.size(); }
  
  /// Returns a iterable sequence of all configured PMT thresholds.
  auto ADCthresholds() const { return util::get_elements<0U>(fADCthresholds); }
  
  /// Returns the ADC threshold value for PMT threshold index `iThr`.
  /// @throws std::out_of_range if `iThr` is not a valid threshold index
  /// @see `nADCthresholds()`
  icarus::trigger::ADCCounts_t ADCthreshold(std::size_t iThr) const
    { return next(fADCthresholds.begin(), iThr)->first; }
  
  /// Returns the detector geometry service provider.
  geo::GeometryCore const& geometry() const { return fGeom; }
  
  /// Returns the detector timings helper.
  detinfo::DetectorTimings const& detTimings() const { return fDetTimings; }
  
  /// Returns the resolution of trigger timing clock [ns]
  nanoseconds triggerTimeResolution() const { return fTriggerTimeResolution; }
  
  /// Returns the first and after-last tick of the beam gate in optical scale.
  std::pair
    <detinfo::timescales::optical_tick, detinfo::timescales::optical_tick>
    const&
    beamGateTickRange() const { return fBeamGateOpt; }
  
  // --- END Helper interface --------------------------------------------------
  
  
  using PlotSandbox = icarus::trigger::PlotSandbox; ///< Import type.
  
  /// List of references to plot sandboxes.
  using PlotSandboxRefs_t
    = std::vector<std::reference_wrapper<PlotSandbox const>>;
  
  
  //----------------------------------------------------------------------------
  class PlotCategory {
    
      public:
    
    /// Type of test function.
    using QualifyFunc_t = std::function<bool(EventInfo_t const&)>;
    
    
    /// Constructor from category name and test function.
    PlotCategory(
      std::string name, std::string descr = {},
      QualifyFunc_t&& test = [](EventInfo_t const&){ return true; }
      )
      : fName(std::move(name)), fDescr(std::move(descr)), fTest(std::move(test))
      {}
    
    /// Returns the name of the category.
    std::string const& name() const { return fName; }
    
    /// Returns the description of the category.
    std::string const& description() const { return fDescr; }
    
    /// Returns whether the event belong to this category.
    bool test(EventInfo_t const& info) const { return fTest(info); }
    
    operator std::string() const { return name(); }
    bool operator() (EventInfo_t const& info) const { return test(info); }
    
      private:
    
    std::string fName;
    std::string fDescr;
    QualifyFunc_t fTest;
    
  }; // class PlotCategory
  using PlotCategories_t = std::vector<PlotCategory>;

  
  //----------------------------------------------------------------------------

  class HistGetter { // helper, since this seems "popular"
    PlotSandbox const& plots;
    
      public:
    HistGetter(PlotSandbox const& plots): plots(plots) {}
    
    PlotSandbox const& box() const { return plots; }
    
    TH1& Hist(std::string const& name) const { return plots.demand<TH1>(name); }
    TH2& Hist2D(std::string const& name) const { return plots.demand<TH2>(name); }
    TEfficiency& Eff(std::string const& name) const
      { return plots.demand<TEfficiency>(name); }
    
  }; // class HistGetter
  
  //----------------------------------------------------------------------------
  
  /// Generic description of trigger settings.
  struct SettingsInfo_t {
    
    std::size_t index;       ///< Settings unique index.
    std::string tag;         ///< A tag of the settings (for object names).
    std::string description; ///< A description of the settings (for plots).
    
    SettingsInfo_t() = default;
    SettingsInfo_t(std::size_t index, std::string tag, std::string descr)
      : index(index), tag(tag), description(descr) {}
    
  }; // SettingsInfo_t
  
  /// Type of trigger gate extracted from the input event.
  using InputTriggerGate_t = icarus::trigger::MultiChannelOpticalTriggerGate;
  
  /// A list of trigger gates from input.
  using TriggerGates_t = std::vector<InputTriggerGate_t>;

  /// Type of lists of gates, one list per cryostat (outer index: cryostat no).
  using TriggerGatesPerCryostat_t = std::vector<TriggerGates_t>;
  
  /// Type of gate data without channel information.
  using TriggerGateData_t = InputTriggerGate_t::GateData_t;
  
  
  // --- BEGIN Customization interface -----------------------------------------
  
  /// Initializes all the plot sets, one per PMT threshold.
  virtual void initializePlots
    (PlotCategories_t categories, std::vector<SettingsInfo_t> const& settings);

  /// Initializes sets of default plots, one per PMT threshold.
  void initializePlots(std::vector<SettingsInfo_t> const& settings)
    { initializePlots(DefaultPlotCategories, settings); }

  /// Initializes full set of plots for (ADC threshold + category) into `plots`.
  virtual void initializePlotSet
    (PlotSandbox& plots, std::vector<SettingsInfo_t> const& settings) const;

  /// Initializes set of plots per complete trigger definition into `plots`.
  virtual void initializeEfficiencyPerTriggerPlots(PlotSandbox& plots) const;
  
  /// Initializes a single, trigger-independent plot set into `plots`.
  virtual void initializeEventPlots(PlotSandbox& plots) const;
  
  /// Returns whether an event with the specified information should be included
  /// in the plots at all (it's a filter).
  virtual bool shouldPlotEvent(EventInfo_t const& eventInfo) const;
  
  /**
   * @brief Simulates all triggers for a trigger setting and plots the results.
   * @param thresholdIndex the index of the PMT threshold of input primitives
   * @param gates the trigger primitives used to simulate the trigger response
   * @param eventInfo general information about the event being simulated
   * @param selectedPlots list of boxes containing plots to be filled
   * 
   * This pure virtual function is the core of the customization of this class.
   * The `simulateAndPlot()` method is expected to:
   * 
   * 1. combine the trigger primitives
   * 2. apply the beam gate
   * 3. generate the trigger response
   * 4. fill all plots
   * 
   * **for each trigger setting**.
   * Note that this helper class has no knowledge of those settings: how many
   * they are, what they mean, which are their parameters. It is
   * `simulateAndPlot()` task to properly process all of them _at once_, i.e.
   * within a single call.
   * 
   * Big deal! The second and fourth step (application of beam gate and plots
   * filling) have some helpers in this same class:
   * 
   * * `applyBeamGateToAll()` can easily apply the beam gate at the right time
   *   (the right time is chosen by the implementation of `simulateAndPlot()`);
   * * filling of standard plots is performed with `fillAllEfficiencyPlots()`,
   *   `fillEfficiencyPlots()` and `fillEventPlots()`, which can be used as such
   *   or be customized for more plots.
   * 
   * Combination of primitive can be helped with
   * `icarus::trigger::TriggerGateData`
   * 
   * All plots must have already been initialized via `initializePlots()` and
   * its helpers.
   * 
   * The method is required to perform the simulation of the trigger primitives
   * specified in input. These primitives are provided in  `gates` as many
   * collections, one per cryostat. While this interface splits the primitives
   * by cryostat, the trigger simulation gets to choose how to handle that
   * information, e.g. performing an independent simulation per cryostat and
   * then combining them all, or just merging back the two collections and
   * treating the two cryostats as a single entity.
   * 
   * The input trigger primitives are extracted from PMT waveforms with a
   * certain PMT threshold, which is documented by its index `thresholdIndex`;
   * the threshold value can be obtained from `ADCthreshold(thresholdIndex)`,
   * and in the future more information may become available as well.
   * 
   * Plots will need information from the event being simulated: that
   * information is precooked, extracted from the event and stored in
   * `eventInfo`. If more information is needed, `EventInfo_t` needs to be
   * updated in this helper, as no way is provided to perform that customization
   * in the derived class.
   * 
   * It is expected that all PMT thresholds and all trigger settings
   * 
   * settings identified
   * by the specified `settings` index. This helper class does not know the
   * definition of any trigger, and it is expected that settings are tracked by
   * the derived classes and each is associated to an index. 
   * 
   */
  virtual void simulateAndPlot(
    std::size_t const thresholdIndex,
    TriggerGatesPerCryostat_t const& gates,
    EventInfo_t const& eventInfo,
    PlotSandboxRefs_t const& selectedPlots
    ) const = 0;
  
  
  /// Fills the plots (`initializeEventPlots()`) with info from `eventInfo`.
  virtual void fillEventPlots
    (EventInfo_t const& eventInfo, PlotSandbox const& plots) const;

  /// Fills the plots (`initializeEfficiencyPerTriggerPlots()`) with info from
  /// `eventInfo` and `triggerInfo`.
  virtual void fillEfficiencyPlots(
    EventInfo_t const& eventInfo,
    TriggerInfo_t const& triggerInfo,
    PlotSandbox const& plots
    ) const;

  /// Fills the plots as `fillEfficiencyPlots()` and also as `fillEventPlots()`
  /// for the proper category: triggered or not triggered.
  virtual void fillAllEfficiencyPlots(
    EventInfo_t const& eventInfo,
    TriggerInfo_t const& triggerInfo,
    PlotSandbox const& plots
    ) const;


  /**
   * @brief Extracts from `event` the relevant information on the physics event.
   *
   * This returns a `EventInfo_t` object filled as completely as possible.
   * The result is used in the context of
   * @ref TriggerEfficiencyPlots_AddingCategory "assigning the `event` to categories".
   *
   * This method can read any data product in the event, and has access to the
   * module configuration. If information is needed from a data product that
   * is not read yet, the following actions are needed:
   *
   * 1. decide the input tag of the data product(s) to be read; this is usually
   *    delegated to the user via module configuration (see `Config` structure);
   * 2. the required input tags need to be stored into `TriggerEfficiencyPlots`
   *    object;
   * 3. _art_ needs to be informed of the intention of reading the data products
   *    via `consumes()` or equivalent calls, in `TriggerEfficiencyPlots`
   *    constructor;
   * 4. the data product can be read with the usual means in this method.
   *
   * For an example of this procedure, see how the configuration parameter
   * `DetectorParticleTag` is treated: read by `Config::DetectorParticleTag`,
   * stored in `fDetectorParticleTag`, declared as "consumed" in
   * `TriggerEfficiencyPlots` constructor. Likewise a known list of data
   * products can be used: see `GeneratorTags` parameter, read by the sequence
   * `Config::GeneratorTags`, stored in the class field `fGeneratorTags`,
   * declared as "consumed" in a loop in `TriggerEfficiencyPlots` constructor,
   * and read in `extractEventInfo()` to be used within it.
   *
   */
  virtual EventInfo_t extractEventInfo(art::Event const& event) const;
  
  // --- END Customization interface -------------------------------------------
  
  
  // --- BEGIN Additional helper utilities -------------------------------------
  
  /// Applies the beam gate coincidence to the specified trigger primitive.
  template <typename GateObject>
  GateObject applyBeamGate(GateObject gate) const;

  /// Applies the beam gate coincidence to all trigger primitives in `gates`.
  template <typename GateColl, typename GateObj = typename GateColl::value_type>
  std::vector<GateObj> applyBeamGateToAll(GateColl const& gates) const;
  
  
  // documentation at definition:
  static PlotCategories_t const DefaultPlotCategories;
  
  
  // --- END Additional helper utilities ---------------------------------------
  
  
    private:
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  std::vector<art::InputTag> fGeneratorTags; ///< Generator data product tags.
  art::InputTag fDetectorParticleTag; ///< Input simulated particles.
  /// Energy deposition data product tags.
  std::vector<art::InputTag> fEnergyDepositTags;
  
  /// ADC thresholds to read, and the input tag connected to their data.
  std::map<icarus::trigger::ADCCounts_t, art::InputTag> fADCthresholds;
  
  /// Duration of the gate during with global optical triggers are accepted.
  microseconds fBeamGateDuration;
  
  nanoseconds fTriggerTimeResolution; ///< Trigger resolution in time.
  
  bool fPlotOnlyActiveVolume; ///< Plot only events in active volume.
  
  /// Message facility stream category for output.
  std::string const fLogCategory;

  std::string fLogEventDetails; ///< Steam where to print event info.
  
  
  /// Plot categories (via `initializePlots()`).
  PlotCategories_t fPlotCategories;
  
  // --- END Configuration variables -------------------------------------------


  // --- BEGIN Service variables -----------------------------------------------

  geo::GeometryCore const& fGeom;
  detinfo::DetectorClocks const& fDetClocks;
  detinfo::DetectorTimings fDetTimings;

  /// ROOT directory where all the plots are written.
  art::TFileDirectory fOutputDir;

  // --- END Service variables -------------------------------------------------

  
  // --- BEGIN Internal variables ----------------------------------------------
  
  /// Gate representing the time we expect light from beam interactions.
  icarus::trigger::OpticalTriggerGate const fBeamGate;
  
  /// Beam gate start and stop tick in optical detector scale.
  std::pair
    <detinfo::timescales::optical_tick, detinfo::timescales::optical_tick>
    const fBeamGateOpt;

  /// Beam gate start and stop time in simulation scale.
  std::pair
    <detinfo::timescales::simulation_time, detinfo::timescales::simulation_time>
    const fBeamGateSim;

  /// ID of cryostat where each optical detector channel (vector index) is.
  std::vector<geo::CryostatID> const fChannelCryostat;
  
  /// All plots, one set per ADC threshold.
  std::vector<PlotSandbox> fThresholdPlots;

  /// Handler of ROOT tree output.
  std::unique_ptr<details::EventIDTree> fIDTree;
  std::unique_ptr<details::PlotInfoTree> fPlotTree; ///< Handler of ROOT tree output.
  std::unique_ptr<details::EventInfoTree> fEventTree; ///< Handler of ROOT tree output.
  
  std::atomic<unsigned int> nEvents { 0U }; ///< Count of seen events.
  std::atomic<unsigned int> nPlottedEvents { 0U }; ///< Count of plotted events.

  // --- END Internal variables ------------------------------------------------


  /// Returns the names of the plot categories event qualifies for.
  std::vector<std::string> selectPlotCategories
    (EventInfo_t const& info, PlotCategories_t const& categories) const;

  /// Returns the TPC `point` is within, `nullptr` if none.
  geo::TPCGeo const* pointInTPC(geo::Point_t const& point) const;
  
  /// Returns in which TPC `point` is within the _active volume_ of;
  /// `nullptr` if none.
  geo::TPCGeo const* pointInActiveTPC(geo::Point_t const& point) const;
  
  /// Reads a set of input gates from the `event`
  /// @return trigger gates, converted into `InputTriggerGate_t`
  TriggerGates_t readTriggerGates
    (art::Event const& event, art::InputTag const& dataTag) const;
  
  /// Moves the data in `gates` in a collection of gates by cryostat.
  TriggerGatesPerCryostat_t splitByCryostat(TriggerGates_t&& gates) const;
  
  
  /// Fills and returns a map of cryostat ID for each optical detector channel.
  static std::vector<geo::CryostatID> makeChannelCryostatMap
    (geo::GeometryCore const& geom);
  
  static std::string thrAndCatName
    (std::string const& boxName, std::string const& category)
    { return boxName + "_" + category; }
  static std::string thrAndCatName
    (PlotSandbox const& box, std::string const& category)
    { return thrAndCatName(box.name(), category); }
  
  
}; // icarus::trigger::TriggerEfficiencyPlotsBase



//==============================================================================
//=== template implementation
//==============================================================================
//--- icarus::trigger::TriggerEfficiencyPlotsBase
//------------------------------------------------------------------------------
template <typename GateObject>
GateObject icarus::trigger::TriggerEfficiencyPlotsBase::applyBeamGate
  (GateObject gate) const
{
  return std::move(gate.Mul(fBeamGate));
} // icarus::trigger::TriggerEfficiencyPlotsBase::applyBeamGate()


//------------------------------------------------------------------------------
template
  <typename GateColl, typename GateObj /* = typename GateColl::value_type */>
std::vector<GateObj>
icarus::trigger::TriggerEfficiencyPlotsBase::applyBeamGateToAll
  (GateColl const& gates) const
{
  
  std::vector<GateObj> res;
  std::transform(
    gates.begin(), gates.end(), std::back_inserter(res),
    [this](GateObj const& gate){ return applyBeamGate(gate); }
    );
  return res;
} // icarus::trigger::TriggerEfficiencyPlotsBase::applyBeamGateToAll()


//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_TRIGGEREFFICIENCYPLOTSBASE_H
