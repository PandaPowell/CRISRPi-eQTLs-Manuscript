#!/bin/bash

# load micromamba into this non-interactive shell
eval "$($MAMBA shell hook -s bash)"
micromamba activate crisprQTL

jupyter nbconvert --to script 08.Intersect_interval_all_crispr.ipynb

"$MAMBA" run -n crisprQTL python "08.Intersect_interval_all_crispr.ipynb"
