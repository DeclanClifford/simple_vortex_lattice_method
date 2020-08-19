# simple_vortex_lattice_method
A proof of concept vortex lattice method for modelling a thin, flat-plate finite-wing. The model is suitable for low speed aerodynamics and accurately models
the effects of changing aspect ratio, sweep angle and dihedral angle

## What is this?
This is a proof of concept vortex lattice method for a flat plate finite wing. 
It is based on an original vortex lattice method by Professor N. Sandham at the University of Southampton. This code goes further than the original introducing wing taper,
sweep and dihedral. The original program and theory this method is based on is given [here](https://github.com/DeclanClifford/simple_vortex_lattice_method/blob/master/theory/vortex-lattice-theory.pdf) 

## Why does this exist?
This program was made with the intention of becoming the foundation of the aerodynamic theory behind a simple flight simulator.
As such it will probably be useless to most people, but is here purely to give evidence of the flight simulator's background.

## How does this work?
See the theory outlined by N. Sandham for a very clear and concise description of the model's theory.

The program, `single_lifting_surface_vlm.m`, outputs the lift coefficient, induced drag coefficient and induced drag factor of a finite wing of user specified geometry.
These outputs
have been chosen to make it easy to verify this model by matching them with the outputs given by N. Sandham in the [theory](https://github.com/DeclanClifford/simple_vortex_lattice_method/blob/master/theory/vortex-lattice-theory.pdf) .

The program `vfil.m` is a dependency of `single_lifting_surface_vlm.m`. A description of this code is also given in the [theory](https://github.com/DeclanClifford/simple_vortex_lattice_method/blob/master/theory/vortex-lattice-theory.pdf) 

It's worth noting the theory states a method to introduce wing twist to the model. I've not added this because the simple flight simulator only uses simple aircraft geomtry,
however, making this addition should be trivial.
