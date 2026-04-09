#!/bin/bash

# TODO: pass path to the script
SIRIL_CLI=siril-cli

WORKDIR=`mktemp -d`

cd $WORKDIR

# download test files
wget https://siril-share-public.s3.rbx.io.cloud.ovh.net/integration/osc/OSC_Preprocessing.ssf
wget https://siril-share-public.s3.rbx.io.cloud.ovh.net/integration/osc/light_1.fit
wget https://siril-share-public.s3.rbx.io.cloud.ovh.net/integration/osc/light_2.fit
wget https://siril-share-public.s3.rbx.io.cloud.ovh.net/integration/osc/light_3.fit
wget https://siril-share-public.s3.rbx.io.cloud.ovh.net/integration/osc/master_dark.fit
wget https://siril-share-public.s3.rbx.io.cloud.ovh.net/integration/osc/master_flat.fit
wget https://siril-share-public.s3.rbx.io.cloud.ovh.net/integration/osc/expected.fit

# simple OSC prepro test
$(SIRIL_CLI) -d `pwd` -s OSC_Preprocessing.ssf 
python3 compare_fits.py expected.fit result_540s.fit

rm -rf $WORKDIR

# TODO: retval should be the return value of python3 compare_fits.py, not of rm
