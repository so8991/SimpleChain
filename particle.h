/*
 
 USC/Viterbi/Computer Science
 "Jello Cube" Assignment 1 starter code
 
 Your name:
 Sunoh Yoo
 
 */
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include "jello.h"
#include <vector>
#include <fstream>

using namespace std;

#ifndef _PARTICLE_H_
#define _PARTICLE_H_
class Particle
{
private:
    void initParticle(double mass);
public:
    //mass of the particle point
    double mass_pnt;
    //position of particle point
    point pnt_position;
    //velocity of particle point
    point pnt_velocity;
    //acceleration of particle point
    point pnt_acceleration;
    
    Particle();
    Particle(double mass);
};
#endif // _PARTICLE_H_
