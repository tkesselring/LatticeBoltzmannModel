# LatticeBoltzmannModel
Implementation of a Lattice Boltzmann Simulation using F#, based on:
[Crude 2D Lattice Boltzmann Demo program](http://www.many-core.group.cam.ac.uk/projects/LBdemo.shtml) by Graham Pullan.

Use build.bat in order to build the solution without Visual studio. 

**Note:** 
The nuget packages are installed using [Paket](https://github.com/fsprojects/Paket).

For the GPU implementation Alea GPU form [QuantAlea](http://quantalea.com/) is used. You can get a 
[Community licence for free](http://quantalea.com/licensing/), it allows you to execute your code on NVidias GeForce GPUs.

Literature:

 - [A Practical Introduction to the Lattice Boltzmann Method](http://www.ndsu.edu/fileadmin/physics.ndsu.edu/Wagner/LBbook.pdf)
 - Collistion terms for different Lattices: [Lattice boltzmann simulations of softmatter systems](http://www.che.ufl.edu/ladd/publications/aps_07.pdf)