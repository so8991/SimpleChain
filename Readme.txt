**** README File **** CSCI 520, Assignment 3
Name: Sunoh Yoo
OS: MacOSX
- In this project, I used particle system to model a simple chain with hard constraint. I used GNU Scientific Library to do this.
- Particle directory contains the project containing source and header files
- Video directory contains a gif video file (movie.gif) and jpg frame images.
<Main Modules> 1. jello.h, jello.cpp
: Those contain main function, contains gl initialization so that it initialize color with lighting information, camera position, etc. I modified the code to declare ParticleSystem object with initialization and its rendering function.
2. particle.h, particle.cpp, particlesystem.h, particlesystem.cpp
: Those contain variables and functions to calculate and render particle system. Particle class contains each particle, and ParticleSystem contains contains the constraint system.
<Acomplished>
- I created Particle and ParticleSystem classes to simulate constraint particle system.
- I set dt = 0.0015, damping coefficient k_damp=0.5, alpha=5, beta = 25/4.
- I set key values to handle forces (left, right, up, down keys and ‘r’:reset, ‘u’: to Euler, ‘k’: to RK4)
<Extra credit>
- I implemented RK4 and Euler integrator
