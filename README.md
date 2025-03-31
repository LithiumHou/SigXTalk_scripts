# SigXTalk_scripts

Scripts for reproducint the results in the article "Dissecting crosstalk induced by cell-cell communication using single-cell transcriptomic data".

## Usage
### Install the dependencies, including:
-- R dependencies:
   Seurat, CellChat, igraph, ranger, tidymodels, tidyfst, parallel, foreach, ggplot2
-- Python dependencies:
   SigXTalk requires a Python module to operate correctly. We recommend that an independent python environment be created to run SigXTalk.
```
conda create -n SigXTalk_py python=3.8
conda activate SigXTalk_py
pip install pandas==2.0.3 scikit-learn==1.3.0 scipy==1.10.1 numpy==1.24.3 argparse==1.4.0
```
SigXTalk could be run on both CUDA and CPU devices. We strongly recommend using CUDA to accelerate the training of neural network using torch:

```
# On Linux or Windows
pip install torch==1.13.1+cu117 torchvision==0.14.1+cu117 torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cu117
# On OSX
pip install torch==1.13.1 torchvision==0.14.1 torchaudio==0.13.1
```
If you do not have a CUDA device, you may use the CPU version of torch. However, it could be quite time-consuming.
```
# On Linux or Windows
pip install torch==1.13.1+cpu torchvision==0.14.1+cpu torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cpu
# On OSX
pip install torch==1.13.1 torchvision==0.14.1 torchaudio==0.13.1
```

Then, the dhg package is required to handle the hypergraph object:
```
pip install dhg
```

### Clone the repository to your device.

### Run any of the demo (e.g., demo_hnscc.R) in the ./scripts directory
