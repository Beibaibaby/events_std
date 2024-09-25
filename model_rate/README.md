# Rate Model for Interleaving Asynchronous and Synchronous Activity

## Overview
This directory contains the MATLAB scripts and C++ code for running the rate model simulations of cortical networks with synaptic plasticity. The model simulates a birth-death process in excitatory and inhibitory populations, exploring the dynamics of synaptic strength modifications over time.

## Files Description
- **markov_2D.cpp**: C++ code for the Markov process simulation, compiled with MATLAB's mex function to be run within MATLAB.
- **run_simulation.m**: MATLAB script that initializes parameters, calls the compiled C++ function, and plots the results.
- **fig6_bifurcation.ode**, **fig7_bifurcation.ode**, **fig8_bifurcation.ode**: Files for generating bifurcation diagrams in XPPAUT, used to analyze the behavior of the system under different parameters.

## Getting Started

### Prerequisites
- **MATLAB**: Ensure you have MATLAB installed on your system, with support for the mex compiler.
- **C++ Compiler**: Ensure your MATLAB installation is configured with a compatible C++ compiler for compiling the mex file.

### Compilation
First, you need to compile the C++ code to be callable from MATLAB:
1. Open MATLAB and navigate to the directory containing `markov_2D.cpp`.
2. Run the following command to compile the code:
   ```matlab
   mex markov_2D.cpp
   ```
   This will generate a `markov_2D.mex*` file (the extension depends on your operating system and MATLAB version).

### Running the Model
To run the simulation:
1. Open `run_simulation.m` with MATLAB.
2. Modify any parameters as necessary to fit your experimental setup.
3. Run `run_simulation.m`. This script sets up initial conditions, calls the compiled Markov process simulation, and plots the activity and plasticity over time.

### Analyzing Bifurcation Diagrams
To generate and analyze bifurcation diagrams:
1. Install XPPAUT on your system, available from [here](http://www.math.pitt.edu/~bard/xpp/xpp.html).
2. Open the `.ode` files (`fig6_bifurcation.ode`, `fig7_bifurcation.ode`, `fig8_bifurcation.ode`) with XPPAUT.
3. Use XPPAUT's analysis tools to explore how the system behavior changes with parameters like synaptic strengths and input drives.

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

---

This README provides all the necessary details to compile, run, and analyze the rate model simulations. It explains how to handle the C++ integration into MATLAB and the use of external tools for further analysis. If you need any further adjustments or additional information, let me know!
