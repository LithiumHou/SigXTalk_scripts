#!/bin/bash
source /home/jiawen/anaconda3/etc/profile.d/conda.sh
echo "Start simulation!"
cd /home/jiawen/myMLnet/benchmark/Beeline/outputs
echo "Round,seed,AUROC,AUPRC" > results_randomseed.csv

# conda activate myMLnet
# Rscript /home/jiawen/myMLnet/benchmark/codes/Simulation.R 10 40 150 500

echo "Simulation finished, now start pathway prediction using SigXTalk..."
conda activate SigXTalk
for cround in $(seq 1 10); do


    cd
    cd myMLnet/benchmark/codes
    python3 /home/jiawen/myMLnet/benchmark/codes/bench_seed.py --round ${cround}


done