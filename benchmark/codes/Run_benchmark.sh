#!/bin/bash
source /home/jiawen/anaconda3/etc/profile.d/conda.sh
echo "Start simulation!"
cd /home/jiawen/myMLnet/benchmark/Beeline/outputs


for cround in $(seq 10 12); do
    find home/jiawen/myMLnet/benchmark/Beeline/inputs/example/bench_rectf -mindepth 1 -type d -exec rm -rf {} +
    find home/jiawen/myMLnet/benchmark/Beeline/inputs/example/bench_tftg -mindepth 1 -type d -exec rm -rf {} +
    rm -rf /home/jiawen/myMLnet/benchmark/Beeline/outputs/example/bench_rectf/GENIE3
    rm -rf /home/jiawen/myMLnet/benchmark/Beeline/outputs/example/bench_tftg/GENIE3
    
    conda activate myMLnet
    Rscript /home/jiawen/myMLnet/benchmark/codes/Simulation.R 10 40 150 500

    echo "Simulation finished, now start pathway prediction using SigXTalk..."
    conda activate SigXTalk

    cd
    cd myMLnet/benchmark/codes
    python /home/jiawen/myMLnet/benchmark/codes/bench_test.py --round ${cround}
    echo "Start prediction using GENELINK..."
    python /home/jiawen/myMLnet/benchmark/codes/bench_genelink.py --round ${cround}


    cd
    cd myMLnet/benchmark/Beeline
    conda activate BEELINE

    echo "Start Rec-TF prediction using other methods..."
    python BLRunner.py --config config-files/config_rectf.yaml
    echo "Start TF-TG prediction using other methods..."
    python BLRunner.py --config config-files/config_tftg.yaml

    cd
    python /home/jiawen/myMLnet/benchmark/codes/bench_other.py ${cround}


done





