# Interleaving Asynchronous and Synchronous Activity in Balanced Cortical Networks with Short Term Synaptic Depression

## Description
This GitHub repository contains the computational models and simulations used in our study, which investigates the dynamics of asynchronous and synchronous activity in cortical networks influenced by short term synaptic depression. Our findings are could be found in:

Dunworth, J. B. *, Xu, Y. *, Graupner, M., Ermentrout, B., Reyes, A. D., & Doiron, B. (2024). Interleaving asynchronous and synchronous activity in balanced cortical networks with short term synaptic depression. TBD

The repository is structured into directories for the rate model (MATLAB) and the spiking model (Julia), featuring scripts that generate traces, run simulations, and plot results.

## Repository Structure
- **model_rate/**: MATLAB scripts for the rate-based model simulations.
  - `markov_2D.cpp`: C++ code compiled with MATLAB’s mex to run Markov processes.
  - `markov_example_generate_traces.m`: MATLAB script that calls the C++ simulation and plots example traces.
  - `fig6_bifurcation.ode`, `fig7_bifurcation.ode`, `fig8_bifurcation.ode`: Files for generating bifurcation diagrams using XPPAUT.
- **model_spiking/**: Julia scripts for the spiking model simulations.
  - `LIF_1.0.0_compare.jl`
  - `run.jl`
  - `run.sbatch`
  - `sim.jl`
  - `hyperparameter/`
- **.DS_Store**: System file for Mac OS folder settings (can be ignored).

## Getting Started

### Prerequisites
- **Julia**: Required for `.jl` scripts. Download from [Julia's official website](https://julialang.org/downloads/).
- **MATLAB or GNU Octave**: Needed for `.m` scripts. MATLAB is also required for compiling and running the `.cpp` file using mex.
- **XPPAUT**: Necessary for running `.ode` files to generate bifurcation diagrams.

### Installation
Clone the repository:
```bash
git clone https://github.com/yourusername/your-repository-name.git
```

### Running the Rate Model (MATLAB)
1. **Compile the C++ code** for MATLAB:
   ```bash
   mex markov_2D.cpp
   ```
2. **Navigate** to the `model_rate` directory.
3. **Execute the MATLAB scripts** to run simulations and plot data:
   ```matlab
   run('markov_example_generate_traces.m')
   ```

### Generating Bifurcation Diagrams (XPPAUT)
- Run the `.ode` files in XPPAUT to generate the bifurcation diagrams for figures 6, 7, and 8.

### Running the Spiking Model (Julia)
1. **Navigate** to the `model_spiking` directory.
2. **Install Julia packages** if not already installed:
   ```julia
   using Pkg
   Pkg.add(["Statistics", "Plots", "Dates", "JSON3", "Random", "JLD2", "Distributed", "SharedArrays", "Measures", "Distributions"])
   ```
3. **Execute the main simulation** script:
   ```bash
   julia run.jl
   ```
   For HPC environments, modify `run.sbatch` as necessary and submit it:
   ```bash
   sbatch run.sbatch
   ```

## Contributing
Contributions are welcome. Please fork the repository and submit pull requests, or open an issue to discuss potential changes.

## License
This project is licensed under the MIT License - see the `LICENSE` file for details.

## Citation
If this model or code aids your research, please cite:
```bibtex
@article{dunworth2024interleaving,
  title={Interleaving asynchronous and synchronous activity in balanced cortical networks with short term synaptic depression},
  author={Dunworth, Jeffrey B.* and Xu, Yunlong* and Graupner, Michael and Ermentrout, Bard and Reyes, Alex D. and Doiron, Brent},
  journal={TBD},
  year={2024},
}
```
