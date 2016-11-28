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
#include "particle.h"

using namespace std;

#ifndef _PARTICLESYSTEM_H_
#define _PARTICLESYSTEM_H_
class ParticleSystem
{
public:
    int no_particle;    //number of particles
    vector<Particle> particle_lst;
    char integrator;
    double dt;
    ofstream error_value;
    int frame;
    //ui
    double left_ext_force;
    double right_ext_force;
    double up_ext_force;
    double down_ext_force;
    //row and col for matrix
    int row;
    int col;
    //len of the link between particles
    double len_link;
    //gravity
    point gravity;
    //constant damping k
    double damp_k;
    //constructor
    ParticleSystem(int num);
    ~ParticleSystem();
    //solve the equation to get acceleration
    void computeAcceleration();
    //particle system renderer
    void glRender(float s_radius, int s_subdivision, float c_radius, int c_subdivision);
    //print out error
    void printError();
    
    
private:
    double link_c_q[4];
    double ring_c_q[2];
    double link_cDot_q[4];
    double ring_cDot_q[2];
    
    //== left matrix ==//
    double * mass_matrix; //Matrix
    double * gradC_matrix; //dc / dq
    double * gradCDot_matrix; //dc' / dq
    double * trans_gradC_matrix; //(dc/dq)T
    double * left_matrix;
    
    //== right matrix ==//
    double * ext_force_matrix; //f(t)
    double *qDot_matrix; //q'
    double * gradCq_DotDot_matrix1; //-(dc'/dq)q'
    double * gradCq_DotDot_matrix2; //(dc/dq)
    double * gradCq_DotDot_matrix3; //c
    double alpha;
    double beta;
    double * right_matrix;
    
    gsl_vector * x;
    
    //== ring variables ==//
    int no_major; //number of major
    int no_minor; //number of minor
    point ring_pos;
    float in_radius; //minor radius
    float out_radius; //major radius
    
    //== compute .. ==//
    double linkConstraint(int idx);
    double ringConstraint(int idx);
    void linkGradient(int idx);
    void ringGradient(int idx);
    void computeLinkCDot(int idx);
    void computeRingCDot(int idx);
    
    void massMatrix();
    void gradCMatrix(); //gradient c matrix
    void gradCDotMatrix(); //gradient c' matrix
    
    void transGradCMatrix(); // transposed gradient c matrix
    void combineToLeft();  //combine matrices to left matrix
    void computeRight();  //compute right matrix
    void extForceMatrix(); //external force matrix
    void matrixMultiple(double * top_matrix, double * new_matrix, double * out_matrix);
    
    //render..
    void glDrawSphere(float radius, int subdivision);
    void glDrawCylin(float radius, int subdivision);
    void glDrawRing();
    
};
#endif
