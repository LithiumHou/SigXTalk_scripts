#!/bin/bash
cd
source /home/jiawen/anaconda3/etc/profile.d/conda.sh
echo "Start simulation!"
cd
cd /home/jiawen/myMLnet/benchmark/Beeline/outputs
echo "Round,LR,batch,AUROC,AUPRC" > results_lr.csv

# test for lr and batch size
lr_array=(0.0001 0.0005 0.001 0.005 0.01 0.02 0.05 0.1)
bs_array=(16 32 48 64 96 128 192 256)

for cround in $(seq 1 5); do
    conda activate myMLnet
    Rscript /home/jiawen/myMLnet/benchmark/codes/Simulation.R 10 40 150 500
    echo "Simulation finished, now start pathway prediction using SigXTalk..."

    for lrt in "${lr_array[@]}" 
    do
        conda activate SigXTalk
        cd
        cd /home/jiawen/myMLnet/benchmark/Beeline/codes
        python3 /home/jiawen/myMLnet/benchmark/codes/bench_lr.py --round "$cround" --lr "$lrt" --batch_size 64
    done

    for bst in "${bs_array[@]}" 
    do
        conda activate SigXTalk
        cd
        cd /home/jiawen/myMLnet/benchmark/Beeline/codes
        python3 /home/jiawen/myMLnet/benchmark/codes/bench_lr.py --round "$cround" --lr 0.01 --batch_size "$bst"
    done
done

