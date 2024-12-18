#!/bin/bash
source /home/jiawen/anaconda3/etc/profile.d/conda.sh
cd
echo "Start simulation!"
cd /home/jiawen/myMLnet/benchmark/Beeline/outputs
echo "Round,ts,AUROC,AUPRC" > results_trainsize.csv

# test for cells

for cround in $(seq 1 10); do

    # source /home/jiawen/anaconda3/etc/profile.d/conda.sh
    conda activate myMLnet
    cd
    Rscript /home/jiawen/myMLnet/benchmark/codes/Simulation.R 10 40 150 500
    echo "Simulation finished, now start pathway prediction using SigXTalk..."
    for round2 in $(seq 1 8); do
        ts=$(echo "0.1 * $round2 + 0.1" | bc)
        conda activate SigXTalk
        cd
        cd myMLnet/benchmark/codes
        python /home/jiawen/myMLnet/benchmark/codes/bench_genes.py --round "$cround" --filen "/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_trainsize.csv" --train_size "$ts"
    done
done