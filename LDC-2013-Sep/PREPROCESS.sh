#!/bin/bash

# Preprocess the data
./cmd.amr-release-proxy.train
./cmd.amr-release-proxy.train.snt
./cmd.amr-release-proxy.train.snt.tok
./cmd.amr-release-proxy.train.tok

# Run the aligner
./cmd.amr-release-proxy.train.aligned
# Remove opN
./cmd.amr-release-proxy.train.aligned.no_opN

# Tag with IllinoisNer
#ssh 192.168.1.3 (supercomputer)
#cd ~/tools/IllinoisNerExtended-v2.3 (remember this is a modified version that outputs CONLL format)
#./tagAMR.train.sh
#./tagAMR.dev.sh
#./tagAMR.test.sh
