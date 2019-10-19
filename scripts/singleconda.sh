#!/usr/bin/env bash
# properties = {properties}
conda activate hla-kir-imp
cyvcf2_installed=$(python -c "import cyvcf2" && echo $?)
matplotlib_installed=$(python -c "import matplotlib" && echo $?)

if
{exec_job}
