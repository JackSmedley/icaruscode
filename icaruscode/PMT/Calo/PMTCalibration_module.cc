////////////////////////////////////////////////////////////////////////
// Class:       PMTCalibration
// Plugin Type: analyzer (art v3_05_00)
// File:        PMTCalibration_module.cc
//
// Generated at Mon Sep 21 15:21:37 2020 by Andrea Scarpelli
//
//  mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////

#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "artdaq-core/Data/Fragment.hh" // Fragment
#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh" // Fragment
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h" // Channel map

#include "canvas/Utilities/Exception.h"

#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"


#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "lardataobj/RawData/OpDetWaveform.h"

#include "icaruscode/PMT/Calo/CaloTools/Waveform.h"
#include "icaruscode/PMT/Calo/CaloTools/Utils.h"

#include "TTree.h"

namespace pmtcalo {
  class PMTCalibration;
}


class pmtcalo::PMTCalibration : public art::EDAnalyzer {

public:

  explicit PMTCalibration(fhicl::ParameterSet const& pset);

  PMTCalibration(PMTCalibration const&) = delete;
  PMTCalibration(PMTCalibration&&) = delete;
  PMTCalibration& operator=(PMTCalibration const&) = delete;
  PMTCalibration& operator=(PMTCalibration&&) = delete;

  //virtual void respondToOpenInputFile(const art::FileBlock& fb) override;
  virtual void beginJob() override;
  virtual void beginSubRun(const art::SubRun &sr) override;

  void getTriggerTime( std::vector<uint16_t> waveform, double & tch16 );

  void analyze(art::Event const& event) override;

  void clean();

private:

  art::InputTag m_fragment_label;
  art::InputTag m_data_label;
  std::string m_channel_dbase;
  bool m_filter_noise;
  fhicl::ParameterSet m_waveform_config;
  size_t m_reference_digitizer_index;

  int m_run;
  int m_event;

  TTree *m_geo_ttree;
  TTree *m_pulse_ttree;

  const icarusDB::IICARUSChannelMap* fChannelMap;

  std::vector<float> *m_channel_id = NULL;
  std::vector<float> *m_baseline = NULL;
  std::vector<float> *m_rms = NULL;
  std::vector<float> *m_peak_time = NULL;
  std::vector<float> *m_amplitude = NULL;
  std::vector<float> *m_integral = NULL;
  std::vector<float> *m_total_charge = NULL;
  std::vector<int>   *m_crate_id = NULL;
  std::vector<float> *m_fragment_id = NULL;
  std::vector<float> *m_fragment_timestamp = NULL;
  std::vector<float> *m_fragment_nseconds = NULL;
  std::vector<float> *m_crate_starttime = NULL;

  // fitted quantities
  std::vector<float> *m_fit_start_time = NULL;
  std::vector<float> *m_error_start_time = NULL;
  std::vector<float> *m_fit_sigma = NULL;
  std::vector<float> *m_error_sigma = NULL;
  std::vector<float> *m_fit_mu = NULL;
  std::vector<float> *m_error_mu = NULL;
  std::vector<float> *m_fit_amplitude = NULL;
  std::vector<float> *m_error_amplitude = NULL;
  std::vector<float> *m_chi2 = NULL;
  std::vector<float> *m_ndf = NULL;
  std::vector<float> *m_fitstatus = NULL; // O:good, >0: bad,  < 0: not working

  art::ServiceHandle<art::TFileService> tfs;

  Waveform *myWaveformAna;

};


//------------------------------------------------------------------------------


pmtcalo::PMTCalibration::PMTCalibration(fhicl::ParameterSet const& pset)
  : art::EDAnalyzer(pset)  // ,
{ 

   m_fragment_label = pset.get<art::InputTag>("FragmentLabel", "daq:CAENV1730");
   m_data_label = pset.get<art::InputTag>("InputModule", "daqPMT");
   m_filter_noise = pset.get<bool>("FilterNoise", false);
   m_waveform_config = pset.get<fhicl::ParameterSet>("WaveformAnalysis");
   m_reference_digitizer_index = pset.get<size_t>("ReferenceDigitizer");

   myWaveformAna = new Waveform(m_waveform_config);

   // Configure the channel mapping services
   fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();

   if( m_reference_digitizer_index < 0 || m_reference_digitizer_index >= 3 ){
    std::cout << "Wrong digitizer index! digitizer index refers to the position of a digitzer in a VME crate. Range is from [0-2]" << std::endl;
    throw;
   }

}


//------------------------------------------------------------------------------


void pmtcalo::PMTCalibration::beginJob()
{

  //For direct light calibration and timing
  m_pulse_ttree = tfs->make<TTree>("pulsetree","tree with laser pulse characterization");

  m_pulse_ttree->Branch("run", &m_run, "run/I" );
  //m_pulse_ttree->Branch("subrun", &m_subrun, "subrun/I" );
  m_pulse_ttree->Branch("event", &m_event, "event/I" );
  m_pulse_ttree->Branch("channel_id", &m_channel_id );
  m_pulse_ttree->Branch("baseline", &m_baseline );
  m_pulse_ttree->Branch("rms", &m_rms );
  m_pulse_ttree->Branch("peak_time", &m_peak_time );
  m_pulse_ttree->Branch("amplitude", &m_amplitude );
  m_pulse_ttree->Branch("integral", &m_integral );
  m_pulse_ttree->Branch("total_charge", &m_total_charge );
  m_pulse_ttree->Branch("m_crate_id", &m_crate_id );
  m_pulse_ttree->Branch("m_fragment_id", &m_fragment_id );
  m_pulse_ttree->Branch("m_fragment_timestamp", &m_fragment_timestamp );
  m_pulse_ttree->Branch("m_fragment_nseconds", &m_fragment_nseconds );
  m_pulse_ttree->Branch("m_crate_starttime", &m_crate_starttime );

  m_pulse_ttree->Branch("fit_start_time", &m_fit_start_time );
  m_pulse_ttree->Branch("error_start_time", &m_error_start_time );
  m_pulse_ttree->Branch("fit_sigma", &m_fit_sigma);
  m_pulse_ttree->Branch("error_sigma", &m_error_sigma);
  m_pulse_ttree->Branch("fit_mu", &m_fit_mu);
  m_pulse_ttree->Branch("error_mu", &m_error_mu);
  m_pulse_ttree->Branch("fit_amplitude", &m_fit_amplitude);
  m_pulse_ttree->Branch("error_amplitude", &m_error_amplitude);
  m_pulse_ttree->Branch("chi2", &m_chi2);
  m_pulse_ttree->Branch("ndf", &m_ndf);
  m_pulse_ttree->Branch("fitstatus", &m_fitstatus);

  m_geo_ttree = tfs->make<TTree>("geotree","tree with detector geo info");

  std::vector<double> pmtX, pmtY, pmtZ;
  std::vector<double> minX, minY, minZ;
  std::vector<double> maxX, maxY, maxZ;
  auto const geop = lar::providerFrom<geo::Geometry>();
  double PMTxyz[3];
  for(size_t opch=0; opch<geop->NOpChannels(); ++opch) {
    geop->OpDetGeoFromOpChannel(opch).GetCenter(PMTxyz);
    pmtX.push_back(PMTxyz[0]);
    pmtY.push_back(PMTxyz[1]);
    pmtZ.push_back(PMTxyz[2]);
  }
  for(auto iter=geop->begin_TPC(); iter!=geop->end_TPC(); ++iter) {
    auto const& tpc = (*iter);
    minX.push_back(tpc.BoundingBox().MinX());
    minY.push_back(tpc.BoundingBox().MinY());
    minZ.push_back(tpc.BoundingBox().MinZ());
    maxX.push_back(tpc.BoundingBox().MaxX());
    maxY.push_back(tpc.BoundingBox().MaxY());
    maxZ.push_back(tpc.BoundingBox().MaxZ());
  }
  m_geo_ttree->Branch("pmtX",&pmtX);
  m_geo_ttree->Branch("pmtY",&pmtY);
  m_geo_ttree->Branch("pmtZ",&pmtZ);
  m_geo_ttree->Branch("minX",&minX);
  m_geo_ttree->Branch("minY",&minY);
  m_geo_ttree->Branch("minZ",&minZ);
  m_geo_ttree->Branch("maxX",&maxX);
  m_geo_ttree->Branch("maxY",&maxY);
  m_geo_ttree->Branch("maxZ",&maxZ);
  
  m_geo_ttree->Fill();

}

//------------------------------------------------------------------------------


 void pmtcalo::PMTCalibration::beginSubRun(const art::SubRun &sr)
 {

   m_run = sr.id().run();
   //m_subrun = 0;//pmtcalo::fileProgNumber(m_filename);

  } // end beginSubRun


//-----------------------------------------------------------------------------


void pmtcalo::PMTCalibration::getTriggerTime( std::vector<uint16_t> waveform, double & tch16 )
{

    //Get the start and end point
    // Find the max val and the maxbin
    // Start point, bin before the one over the threshold (5% max)
    // End point first bin on the plateau

    double max=0;
    double binmax=0;

    for( size_t bin=0; bin<waveform.size(); bin++ ){

      double value = waveform[bin];
      if( value > max ){
        max = value;
        binmax=bin;
      }
    }

    // Now we found the beam start
    double startval=0;
    double binstart=0;

    for( size_t bin=0; bin<waveform.size(); bin++ ){

      double value = waveform[bin+1]-waveform[bin];
      double delta = max - waveform[bin+1];


      if( value>0.05*delta ){
        startval = waveform[bin];
        binstart=bin;
        break;
      }
    }

    // Now we found the beam end
    double endval=0;
    double binend=0;

    for( size_t bin=binstart; bin<binmax+1; bin++ ){

      double value = waveform[bin];
      if( max-value < 10  ){
        endval = value;
        binend = bin;
        break;
      }
    }

    //std::cout << binmax << " " << max << " " << binmax*2+1 << std::endl;
    //std::cout << binstart << " " << startval << " " << binstart*2+1 << std::endl;
    //std::cout << binend << " " << endval << " " << binend*2+1 << std::endl;

    // We now calculate the start time

    double tstart = (double)binstart*2+1;
    double tend = (double)binend*2+1;
    
    int npx = 1000;

    // Linear interpolation y = a + m*x
    double m = (endval-startval)/(tend-tstart);
    double a = endval - m*tend;

    for( int i=0; i<npx; i++ ){

      double delta = endval-startval;

      double t = tstart + (double)i*(tend-tstart)/npx;
      double val =  a+m*t - startval;

      //std::cout << i*(tend-tstart)/npx << " " << t << " " << val << std::endl;

      if( val > 0.05*delta ){
        tch16 = t;
        break;
      }

    }

    return;
}


//-----------------------------------------------------------------------------


void pmtcalo::PMTCalibration::analyze(art::Event const& event)
{

   m_event = event.id().event();


   // Get the association between channelID and fragment

   std::map<raw::Channel_t, size_t > m_frag_map;
   std::map<raw::Channel_t, sbndaq::CAENV1730FragmentMetadata > m_metafrag_map;
   std::map<int, double> m_crate_time;

   std::cout << m_fragment_label << std::endl;

   auto const& daq_handle = event.getValidHandle<artdaq::Fragments>(m_fragment_label);

   if (daq_handle.isValid() && daq_handle->size() > 0) {


      double fragment_start_time=0;


      for (auto const &artdaqFragment: *daq_handle) {

       size_t fragment_id = artdaqFragment.fragmentID();

       sbndaq::CAENV1730Fragment         fragment(artdaqFragment);
       sbndaq::CAENV1730FragmentMetadata metafrag = *fragment.Metadata();
	     sbndaq::CAENV1730Event evt = *fragment.Event();
	     sbndaq::CAENV1730EventHeader header = evt.Header;

       size_t nChannelsPerBoard  = metafrag.nChannels; //fragment.nChannelsPerBoard();
       uint32_t ev_size_quad_bytes         = header.eventSize;
       uint32_t evt_header_size_quad_bytes = sizeof(sbndaq::CAENV1730EventHeader)/sizeof(uint32_t);
       uint32_t data_size_double_bytes     = 2*(ev_size_quad_bytes - evt_header_size_quad_bytes);
       uint32_t nSamplesPerChannel         = data_size_double_bytes/nChannelsPerBoard;

       const uint16_t* data_begin = reinterpret_cast<const uint16_t*>(artdaqFragment.dataBeginBytes() + sizeof(sbndaq::CAENV1730EventHeader));
       const uint16_t* value_ptr  = data_begin;
       uint16_t        value      = 0;
       size_t          ch_offset  = 0;
	
       if (fChannelMap->hasPMTDigitizerID(fragment_id))
       {
          const icarusDB::DigitizerChannelChannelIDPairVec& digitizerChannelVec = fChannelMap->getChannelIDPairVec(fragment_id);
          
          // Get a new start time when we're looking to fragments at the beginning of each vme crate
          if( fragment_id % 3 == m_reference_digitizer_index ){

            // create the waveoform and pre-allocate the space for the sample
            std::vector<uint16_t> wvfm(nSamplesPerChannel);

            ch_offset = 15 * nSamplesPerChannel;

            for(size_t i_t=0; i_t < nSamplesPerChannel; ++i_t)
            {
                value_ptr = data_begin + ch_offset + i_t; /*pointer arithmetic*/
                value     = *(value_ptr);
                wvfm[i_t] = value;
            }

            // Now that the waveform for ch16 is created, we get the starting time
            getTriggerTime( wvfm, fragment_start_time );

            int crate_id = (int)fragment_id / 3;
            m_crate_time[crate_id] = fragment_start_time;

          }
          

          for(const auto& digitizerChannelPair : digitizerChannelVec)
          {
              raw::Channel_t channelID = digitizerChannelPair.second;
              m_metafrag_map[channelID] = metafrag;
              m_frag_map[channelID] = fragment_id;
          }
        } 
      }
    }

   art::Handle< std::vector< raw::OpDetWaveform > > rawHandle;
   event.getByLabel(m_data_label, rawHandle);

   // There is a valid handle per channel
   for( auto const& raw_waveform : (*rawHandle) )
   {

     // Without the correct mapping, we need to set the association with
     // the digitizer id by hand (in future this will be moved to the decoder)

     raw::Channel_t channel_id = raw_waveform.ChannelNumber();

     m_channel_id->push_back( channel_id );
     m_fragment_id->push_back( m_frag_map[channel_id] );
     m_fragment_timestamp->push_back( m_metafrag_map[channel_id].timeStampSec );
     m_fragment_nseconds->push_back( m_metafrag_map[channel_id].timeStampNSec );

     int crate_id = (int)m_frag_map[channel_id] / 3;
     m_crate_id->push_back( crate_id );
     m_crate_starttime->push_back( m_crate_time[crate_id] );

     myWaveformAna->loadData( raw_waveform );
     if( m_filter_noise ){ myWaveformAna->filterNoise(); }

     auto pulse = myWaveformAna->getLaserPulse();

     m_baseline->push_back( myWaveformAna->getBaselineMean() );
     m_rms->push_back( myWaveformAna->getBaselineWidth() );
     m_peak_time->push_back( pulse.time_peak );
     m_amplitude->push_back( pulse.amplitude );
     m_integral->push_back( pulse.integral );
     m_total_charge->push_back( myWaveformAna->getTotalCharge() );

     m_fit_start_time->push_back(pulse.fit_start_time);
     m_error_start_time->push_back(pulse.error_start_time);
     m_fit_sigma->push_back(pulse.fit_sigma);
     m_error_sigma->push_back(pulse.error_sigma);
     m_fit_mu->push_back(pulse.fit_mu);
     m_error_mu->push_back(pulse.error_mu);
     m_fit_amplitude->push_back(pulse.fit_amplitude);
     m_error_amplitude->push_back(pulse.error_amplitude);
     m_chi2->push_back(pulse.chi2);
     m_ndf->push_back(pulse.ndf);
     m_fitstatus->push_back(pulse.fitstatus);

     // Prepare for the next event
     myWaveformAna->clean();

   } // end loop over pmt channels

   m_pulse_ttree->Fill();

   // Cancel the arrays
   clean();


} // end analyze


//-----------------------------------------------------------------------------

void pmtcalo::PMTCalibration::clean(){

  m_channel_id->clear();
  m_baseline->clear();
  m_rms->clear();
  m_peak_time->clear();
  m_amplitude->clear();
  m_integral->clear();
  m_total_charge->clear();
  m_crate_id->clear();
  m_fragment_id->clear();
  m_fragment_timestamp->clear();
  m_fragment_nseconds->clear();
  m_crate_starttime->clear();

  m_fit_start_time->clear();
  m_error_start_time->clear();
  m_fit_sigma->clear();
  m_error_sigma->clear();
  m_fit_mu->clear();
  m_error_mu->clear();
  m_fit_amplitude->clear();
  m_error_amplitude->clear();
  m_chi2->clear();
  m_ndf->clear();
  m_fitstatus->clear();

}


DEFINE_ART_MODULE(pmtcalo::PMTCalibration)