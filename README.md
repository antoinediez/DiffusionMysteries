# The Wonders of Diffusion Phenomena

This repository contains the code used to generate the videos shown in a public lecture given by [Hiroshi Matano](https://researchmap.jp/hmatano?lang=en) (俣野 博) entitled *The wonders of diffusion phenomena* (in Japanese 「拡散現象の不思議」). This lecture was given at the Institute of Statistical Mathematics (統計数理研究所) in Tokyo on 18th November 2023. 

Please refer to [the webpage of the event](https://www.ism.ac.jp/events/2023/meeting1118.html) for more information. A recording of the lecture can be found [here](https://www.ism.ac.jp/events/2023/meeting1118.html) (in Japanese). 

## Structure of this repository

Each `.ipynb` Jupyter notebook file corresponds to one of the experiments described below. These simulations are intended to explain to a general audience how diffusion naturally emerge from discrete random systems with simple rules as in the famous [Brown experiment](https://en.wikipedia.org/wiki/Brownian_motion). The Julia script `RD_solver.jl` is a basic reaction-diffusion solver written in [Julia](https://julialang.org/) which can be used to generate [Turing patterns](https://en.wikipedia.org/wiki/Turing_pattern). 

The technical details and precise implementation are described in the folder `scr`. 

### Hard-sphere motion

The (deterministic) motion of colliding hard-spheres starting from a random initial condition is a simple, yet historical system from which Brownian motion emerges (see for instance [this article](https://doi.org/10.1007/s00222-015-0593-9) for a rigorous but very advanced rigorous mathematical proof of this). 

The two notebooks `hardspheres.ipynb` and `hardspheres_twopopulations.ipynb` show a numerical illustration of this fact. In these simulations, a collision between two hard-spheres is assumed to conserve energy and momentum. In the second case, two systems with a small and a large total energy are shown to illustrate how it affects the behaviour of the systems. 

Example videos are shown below. 

### Lattice motion in 2D

A simplified system consists in a random walk on a lattice: in the 2D case, a particle is allowed to move to one of the four neighbouring positions (up, down, left, right) with equal probability at some random times. After zooming out and accelerating time, the motion of the particle resembles a Brownian motion. 

An example video produced by the notebook `dim2.ipynb` is shown below. 


### Lattice motion in 1D

In a 1D lattice system, particles can move either on the left or on the right (with equal probability). When many (independent) particles are simulated, a distinctive smoothing effect can be observed.

Example videos produced by the notebook `dim1.ipynb` and `dim1_twopopulations.ipynb` are shown below. 

### Turing patterns 

At the continuum level, diffusion can be coupled with "chemical" reactions to form a so-called [reaction-diffusion system](https://en.wikipedia.org/wiki/Reaction%E2%80%93diffusion_system). Since the pioneering work of [Turing](https://en.wikipedia.org/wiki/Turing_pattern), it is known that diffusion can be the source of instabilities leading to patterns such as stripes and spots commonly observed in nature. 

An example of reaction-diffusion system is the so-called Schnakenberg system. Spots and stripes are produced when the diffusion is high enough but disappear when the diffusion vanishes, as shown in the video below produced by the script `RD_solver.jl`. 

## Licence

The videos and the (modest) code can be freely used and modified, with proper citation, under the permissive MIT licence. 

**Authors:** [Antoine Diez](https://antoinediez.gitlab.io/), Hiroshi Matano 