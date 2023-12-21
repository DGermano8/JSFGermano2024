#!/usr/bin/env bash
#
# Run the TIV model on the Baccam data sets for each individual
# separately as described in their manuscript.
#
# This script will create a directory for each patient in the out/
# directory and run the model for each patient.  The output will be
# placed in the patient's directory.
#
for i in {1..6}
do
    echo "Running patient $i"
    if [ ! -d "out/patient-$i" ]; then
	mkdir out/patient-$i
    fi
    python run-inf-cli.py --out out/patient-$i --input_toml config/cli-tiv-jsf.toml --obs_ssv data/patient-$i.ssv --param_plots --state_plot --record_summary
done
