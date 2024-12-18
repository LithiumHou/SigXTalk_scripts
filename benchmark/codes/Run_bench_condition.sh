#!/bin/bash
source /home/jiawen/anaconda3/etc/profile.d/conda.sh

echo "Start simulation!"
cd /home/jiawen/myMLnet/benchmark/Beeline/outputs
echo "Round,Method,AUROC,AUPRC" > results_genes.csv
cd
# test for cells
for ngs in $(seq 1 5); do
    
    for cround in $(seq 1 5); do

        let ng1=5*ngs
        let ng2=20*ngs
        let ng3=75*ngs

        source /home/jiawen/anaconda3/etc/profile.d/conda.sh
        conda activate myMLnet
        Rscript /home/jiawen/myMLnet/benchmark/codes/Simulation.R "$ng1" "$ng2" "$ng3" 500

        echo "Simulation finished, now start pathway prediction using SigXTalk..."
        conda activate SigXTalk
        cd
        cd myMLnet/benchmark/codes
        python /home/jiawen/myMLnet/benchmark/codes/bench_genes.py --round "$cround"
    done
done

#!/bin/bash
source /home/jiawen/anaconda3/etc/profile.d/conda.sh
cd
echo "Start simulation!"
cd /home/jiawen/myMLnet/benchmark/Beeline/outputs
echo "Round,Cells,AUROC,AUPRC" > results_cells.csv


# test for cells
for ncs in $(seq 1 5); do
    let ncells=250*ncs
    for cround in $(seq 1 5); do
  
        conda activate myMLnet
        Rscript /home/jiawen/myMLnet/benchmark/codes/Simulation.R 10 40 150 "$ncells"

        echo "Simulation finished, now start pathway prediction using SigXTalk..."
        conda activate SigXTalk
        cd myMLnet/benchmark/codes
        python /home/jiawen/myMLnet/benchmark/codes/bench_cells.py --round "$cround" 
    done
done

#!/bin/bash
source /home/jiawen/anaconda3/etc/profile.d/conda.sh

echo "Start simulation!"
cd /home/jiawen/myMLnet/benchmark/Beeline/outputs
echo "Round,Method,AUROC,AUPRC" > results_netdense.csv

# test for cells
nds_array=(0.15 0.25 0.35 0.5 0.75 1.0)

for nd in ${nds_array[@]}; do
    
    for cround in $(seq 1 5); do

        cd
        source /home/jiawen/anaconda3/etc/profile.d/conda.sh
        conda activate myMLnet
        cd
        Rscript /home/jiawen/myMLnet/benchmark/codes/Simulation_netP.R "$nd"

        echo "Simulation finished, now start pathway prediction using SigXTalk..."
        conda activate SigXTalk
        cd
        cd myMLnet/benchmark/codes
        python /home/jiawen/myMLnet/benchmark/codes/bench_genes.py --round "$cround" --filen "/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_netdense.csv" --netP "$nd"
    done
done