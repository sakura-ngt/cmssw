#!/bin/bash -ex

#LOCALPATH='/eos/cms/store/group/tsg-phase2/user/jprendi/NERD25/MoreStats/Muons/HLT/'
LOCALPATH='/eos/cms/store/group/tsg/STORM/user/musich/EGamma/HLT/run_403818'
echo "Input source: |${LOCALPATH}|"
LOCALFILES=$(ls -1 ${LOCALPATH} | grep Scouting)
ALL_FILES=""
for f in ${LOCALFILES[@]}; do
    ALL_FILES+="file:${LOCALPATH}/${f},"
done

# Remove the last character
ALL_FILES="${ALL_FILES%?}"
echo "Discovered files: $ALL_FILES"

cmsDriver.py step2 \
	     -s DQM:ScoutingDileptonMonitor \
	     --conditions 160X_dataRun3_HLT_v1 \
	     --data \
	     --datatier DQMIO \
	     -n 100000 \
	     --eventcontent DQMIO \
	     --geometry DB:Extended \
	     --era Run3_2026 \
	     --filein "${ALL_FILES}" \
	     --fileout file:step2.root \
	     --nThreads 24 \
	     --python_filename dqm.py \
	     --no_exec


cat <<@EOF >> dqm.py
process.source.inputCommands = cms.untracked.vstring(
    'keep *',
    'drop CTPPSDiamondDigiedmDetSetVector_hltCTPPSDiamondRawToDigi_TimingDiamond_HLTX',
    'drop CTPPSPixelDigiedmDetSetVector_hltCTPPSPixelDigis__HLTX',
    'drop CTPPSRecord_hltOnlineMetaDataDigis__HLTX',
    'drop DCSRecord_hltOnlineMetaDataDigis__HLTX',
    'drop OnlineLuminosityRecord_hltOnlineMetaDataDigis__HLTX',
    'drop recoBeamSpot_hltOnlineMetaDataDigis__HLTX',
)

from RecoVertex.BeamSpotProducer.BeamSpotOnline_cfi import onlineBeamSpotProducer as _onlineBeamSpotProducer
process.hltOnlineBeamSpot = _onlineBeamSpotProducer.clone()

from PhysicsTools.Scouting.Run3ScoutingElectronBestTrackProducer_cfi import Run3ScoutingElectronBestTrackProducer as _Run3ScoutingElectronBestTrackProducer
process.run3ScoutingElectronBestTrack =  _Run3ScoutingElectronBestTrackProducer.clone()

process.scoutingCollectionMonitor.onlineMetaDataDigis = "hltOnlineMetaDataDigis"
process.scoutingCollectionMonitor.onlyScouting = True
process.scoutingCollectionMonitor.beamSpot = cms.InputTag("hltOnlineBeamSpotFromDB")
process.ScoutingTrackMonitor.beamSpotLabel = cms.InputTag("hltOnlineBeamSpotFromDB")
process.dqmoffline_step = cms.EndPath(process.hltOnlineBeamSpot+process.run3ScoutingElectronBestTrack+process.ScoutingDileptonMonitor)
@EOF

cmsRun dqm.py >& dqm.log

cmsDriver.py step3 -s HARVESTING:@standardDQM \
	     --conditions 160X_dataRun3_HLT_v1 \
	     --data \
	     --geometry DB:Extended \
	     --scenario pp \
	     --filetype DQM \
	     --era Run3_2025 \
	     -n 10000 \
	     --filein file:step2.root \
	     --fileout file:step3.root \
	     --no_exec

cmsRun step3_HARVESTING.py >& harvesting.log
