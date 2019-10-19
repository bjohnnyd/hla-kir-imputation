#!/usr/bin/env bash
#$ -cwd
#$ -V

python -c "import cyvcf2" 2> /dev/null || pip install cyvcf2==0.11.5
python -c "import matplotlib"  2> /dev/null || pip install matplotlib==3.1.1

{exec_job}
