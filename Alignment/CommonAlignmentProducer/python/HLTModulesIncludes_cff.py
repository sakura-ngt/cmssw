import FWCore.ParameterSet.Config as cms

# Run-3 HLT ESProducers, taken from the HLT menu to avoid going out of sync
from HLTrigger.Configuration.HLT_FULL_cff import fragment

hltESPStripCPEfromTrackAngle = fragment.hltESPStripCPEfromTrackAngle
hltESPPixelCPEGeneric        = fragment.hltESPPixelCPEGeneric
hltESPTTRHBWithTrackAngle    = fragment.hltESPTTRHBWithTrackAngle

# ------------------------------------------------------------------------------
# Phase-2: take the equivalents from the 75e33 menu, but keep the Run-3 labels /
# ComponentNames as the interface, so that consumers (TTRHBuilder = 'hltESPTTRHBWithTrackAngle')
# need no era-specific code.
from Configuration.Eras.Modifier_phase2_common_cff import phase2_common


from RecoLocalTracker.Phase2TrackerRecHits.Phase2StripCPEESProducer_cfi import phase2StripCPEESProducer
#from RecoLocalTracker.SiPixelRecHits.PixelCPEGeneric_cfi import PixelCPEGenericESProducer as _hltESPPixelCPEGenericPhase2
#from HLTrigger.Configuration.HLT_75e33.eventsetup.hltESPPhase2StripCPE_cfi import hltESPPhase2StripCPE
#from HLTrigger.Configuration.HLT_75e33.eventsetup.hltESPPixelCPEGeneric_cfi import hltESPPixelCPEGeneric as _hltESPPixelCPEGenericPhase2
from HLTrigger.Configuration.HLT_75e33.eventsetup.hltESPTTRHBuilderWithTrackAngle_cfi import hltESPTTRHBuilderWithTrackAngle as _hltESPTTRHBuilderWithTrackAnglePhase2

# same label, different content
#phase2_common.toReplaceWith(hltESPPixelCPEGeneric, _hltESPPixelCPEGenericPhase2)

# Phase-2 builder content, published under the Run-3 ComponentName
phase2_common.toReplaceWith(hltESPTTRHBWithTrackAngle,
    _hltESPTTRHBuilderWithTrackAnglePhase2.clone(ComponentName = 'hltESPTTRHBWithTrackAngle'))
