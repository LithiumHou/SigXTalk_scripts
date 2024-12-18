#!/bin/bash

# Run the benchmark
echo "Start evaluation!"

# target type: CAF, Endothelial, malignant, myofibroblast, T_cell
cd
# Rscript /home/jiawen/myMLnet/benchmark/codes/benchmark2_analysis.R "CAF" "REACTOME"
# Rscript /home/jiawen/myMLnet/benchmark/codes/benchmark2_analysis.R "Endothelial" "REACTOME"
# Rscript /home/jiawen/myMLnet/benchmark/codes/benchmark2_analysis.R "malignant" "REACTOME"
# Rscript /home/jiawen/myMLnet/benchmark/codes/benchmark2_analysis.R "myofibroblast" "REACTOME"
# Rscript /home/jiawen/myMLnet/benchmark/codes/benchmark2_analysis.R "T_cell" "REACTOME"

Rscript /home/jiawen/myMLnet/benchmark/codes/benchmark2_coexp.R 
Rscript /home/jiawen/myMLnet/benchmark/codes/benchmark2_ppr.R 
