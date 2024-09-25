# Spiking Model for Interleaving Asynchronous and Synchronous Activity

## Overview
This directory contains the Julia scripts and support files for running the spiking model simulations of cortical networks. The model explores the dynamics of neural activity, focusing on the effects of synchronous and asynchronous interactions, influenced by configurable network parameters.

## Files Description
- **LIF_1.0.0_compare.jl**: Contains utility functions used across other scripts, supporting the main simulation processes.
- **run.jl**: Entry point script that sets up the environment and runs the simulations.
- **run.sbatch**: Script for submitting jobs to a high-performance computing cluster.
- **sim.jl**: Core script including network simulation functions.
- **hyperparameter/**:
  - **weights_records**: Stores weight matrices for simulations. A default matrix is provided for consistent experimental setups; however, you can also use custom matrices to explore different configurations.

## Getting Started

### Prerequisites
Ensure Julia is installed on your system. You can download it from [Julia's official website](https://julialang.org/downloads/). Additionally, the following Julia packages are required:

```julia
using Statistics
using Plots  
using Dates  # For generating timestamps
using JSON3  # For JSON handling
using Random
using JLD2
using Distributed
using SharedArrays
using Measures
using Distributions
```

### Installation
First, clone the main project repository:
```bash
git clone https://github.com/yourusername/your-repository-name.git
```
Navigate to the `model_spiking` directory:
```bash
cd your-repository-name/model_spiking
```

Install the required Julia packages by starting Julia and running:
```julia
using Pkg
Pkg.add(["Statistics", "Plots", "Dates", "JSON3", "Random", "JLD2", "Distributed", "SharedArrays", "Measures", "Distributions"])
```

### Running the Model
To run the simulation with default parameters, execute:
```bash
julia run.jl
```

For large-scale simulations or those requiring significant computational resources, modify `run.sbatch` according to your cluster configuration and submit it using:
```bash
sbatch run.sbatch
```

### Hyperparameters and Network Weights
You can modify network parameters and initial conditions by editing the values in the `hyperparameter/weights_records` file or by providing your own.

#### Downloading Custom Weights
If you want to use custom weights for your simulations, you can create a random matrix and store it, or download an example one:
1. Visit the Google Drive link provided: [Download Weights](https://drive.google.com/drive/folders/xxxxxx)
2. Download the `weights_records.jdl2` file.
3. Place the file in the `hyperparameter/` directory.
4. Update the path in `sim.jl` to point to your new weights file, ensuring the format matches that expected by the simulation scripts.

## Contributing
Contributions to this project are welcome. Please fork the repository and submit pull requests with your changes, or open an issue to discuss potential modifications.

## Citation
If this model or the associated code aids your research, please cite our work:
```bibtex
@article{dunworth2024interleaving,
  title={Interleaving asynchronous and synchronous activity in balanced cortical networks with short term synaptic depression},
  author={Dunworth, Jeffrey B. and Xu, Yunlong and Graupner, Michael and Ermentrout, Bard and Reyes, Alex D. and Doiron, Brent},
  journal={TBD},
  year={2024},
}
```
