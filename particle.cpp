/*
 
 USC/Viterbi/Computer Science
 "Jello Cube" Assignment 1 starter code
 
 Your name:
 Sunoh Yoo
 
 */
#include "particle.h"

Particle::Particle()
{
    Particle::initParticle(0.0);
}
Particle::Particle(double mass)
{
    initParticle(mass);
}
void Particle::initParticle(double mass)
{
    mass_pnt = mass;
    pnt_position.x = 0.0;
    pnt_position.y = 0.0;
    pnt_position.z = 0.0;
    
    pnt_velocity.x = 0.0;
    pnt_velocity.y = 0.0;
    pnt_velocity.z = 0.0;
    
    pnt_acceleration.x = 0.0;
    pnt_acceleration.y = 0.0;
    pnt_acceleration.z = 0.0;
}
