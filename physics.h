/*
 
 USC/Viterbi/Computer Science
 "Jello Cube" Assignment 1 starter code
 
 Your name:
 Sunoh Yoo
 
 */

#ifndef _PHYSICS_H_
#define _PHYSICS_H_
#include "particlesystem.h"
//#include "particlesystem.h"
//void computeAcceleration(struct world * jello, struct point a[8][8][8]);

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(ParticleSystem *ps);
void RK4(ParticleSystem *ps);

#endif

