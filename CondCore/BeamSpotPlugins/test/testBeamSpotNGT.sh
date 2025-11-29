#!/bin/bash
# Save current working dir so img can be outputted there later
W_DIR=$(pwd);
source /afs/cern.ch/cms/cmsset_default.sh;
eval `scram run -sh`;
# Go back to original working directory
cd $W_DIR;
# Run get payload data script

mkdir -p $W_DIR/results

getPayloadData.py \
    --plugin pluginBeamSpot_PayloadInspector \
    --plot plot_BeamSpotParametersDiffTwoTags \
    --tag BeamSpotObject_NGT_FullTracking \
    --tagtwo BeamSpotObject_NGT_TwoIterations \
    --input_params "{}" \
    --time_type Run \
    --iovs '{"start_iov": "1", "end_iov": "1"}' \
    --iovstwo '{"start_iov": "1", "end_iov": "1"}' \
    --db Prep \
    --test

mv *.png  $W_DIR/results/beamspot_iters01.png
