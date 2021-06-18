///////////////////////////////////////////////////////////////////////
///
/// \file   IRawDigitFilter.h
///
/// \brief  This provides an interface for tools which are tasked with
///         filtering input waveforms
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IRawDigitFilter_H
#define IRawDigitFilter_H

#include "fhiclcpp/ParameterSet.h"
#include "icaruscode/TPC/SignalProcessing/RawDigitFilter/Algorithms/RawDigitNoiseFilterDefs.h"
#include "TProfile.h"

namespace art
{
    class TFileDirectory;
}

namespace caldata
{
    enum HistogramType : int
    {
        WAVEFORM,
        WAVELESSAVE,
        EROSION,
        DILATION,
        AVERAGE,
        DIFFERENCE,
        OPENING,
        CLOSING,
        DOPENCLOSING,
        LASTELEMENT
    };
    
    using HistogramMap = std::map<int, TProfile*>;

    class IRawDigitFilter
    {
    public:
        virtual ~IRawDigitFilter() noexcept = default;
        
        virtual void   configure(const fhicl::ParameterSet& pset)                      = 0;
        virtual void   initializeHistograms(art::TFileDirectory&)                const = 0;
        virtual size_t plane()                                                   const = 0;
        
        // Find the ROI's
        virtual void FilterWaveform(RawDigitVector&, size_t, size_t, float = 0.) const = 0;
    };
}

#endif
