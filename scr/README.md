# Implementation details 

This folder contains the core scripts used in the simulation experiments. 

## Particle systems 

The simulations using particles systems are written in Python. 


### Hard-sphere system

The hard-sphere system is simulated using the high-performance [SiSyPHE library](https://github.com/antoinediez/Sisyphe) which heavily relies on [PyTorch](https://pytorch.org/) and the [KeOps library](https://www.kernel-operations.io/keops/index.html). 

There is no particular trick to simulate this system. It is discretized in time using fixed (adaptive) time-steps. At each time step, the collision events are solved for all pairs of overlapping spheres. If too many collisions have occured during one time steps or if the system is so crowded that solving a collision creates a new overlapping situation, a smaller time step is chosen. This codes remains relatively fast even for large systems thanks to the high parallelization implementation provided by the [KeOps library](https://www.kernel-operations.io/keops/index.html).

### Random walks 

For the 2D and 1D lattice models, a random walk is implemented as an instance of the class `RandomWalk` which contains most of the methods related to the "mathematics". It can then be included in an instance of the class `Simulation` which contains most of the method to actually plot and render the simulation. The background is plotted and updated as an instance of the class `Lattice`. 

### Reaction-Diffusion systems

The folder `reaction_diffusion_solver.jl` is a basic reaction-diffusion solver which relies heavily on the state-of-the-art [Julia ODE solvers](https://docs.sciml.ai/DiffEqDocs/stable/). 