/*
 
 USC/Viterbi/Computer Science
 "Jello Cube" Assignment 1 starter code
 
 Your name:
 Sunoh Yoo
 
 */
#include "particlesystem.h"

ParticleSystem::ParticleSystem(int num)
{
    no_particle = num;
    len_link = 1.0/(no_particle-1);
    double m = 1.0/(no_particle);
    
    Particle p(m);
    particle_lst.resize(no_particle, p);
    particle_lst[0].mass_pnt = 0.0;
    //first particle --> 0
    
    //== initialize all variables...
    gravity.x = 0.0;
    gravity.y = -1.0;
    gravity.z = 0.0;
    damp_k = 0.5;
    integrator = 'E';
    dt = 0.0015;
    for(int i=0;i<4;i++)
    {
        link_c_q[i] = 0.0;
        link_cDot_q[i] = 0.0;
    }
    for(int i=0;i<2;i++)
    {
        ring_c_q[i] = 0.0;
        ring_cDot_q[i] = 0.0;
    }
    row = no_particle;
    col = 2*(no_particle-1);
    beta = (5.0*5.0)/4.0;
    //ring
    in_radius = 0.0075;
    out_radius = 0.5;
    no_major = 30;
    no_minor = 30;
    ring_pos.x = 0.0;
    ring_pos.y = -0.5;
    ring_pos.z = 0.0;
    frame = 0;
    up_ext_force = 0.0;
    down_ext_force = 0.0;
    left_ext_force = 0.0;
    right_ext_force = 0.0;
    for(int i=0;i<(no_particle-1)/2;i++)
    {
        particle_lst[i].pnt_position.x = 0.0;
        particle_lst[i].pnt_position.y = particle_lst[i].pnt_position.y-len_link*i;
        particle_lst[i].pnt_position.z = 0.0;
    }
    for (int i = (no_particle - 1) / 2; i < no_particle; i++)
    {
        particle_lst[i].pnt_position.x = particle_lst[i].pnt_position.x + len_link * (i-(no_particle-1)/2);
        particle_lst[i].pnt_position.y = -0.5;
        particle_lst[i].pnt_position.z = 0.0;
    }
    // allocate memory to all matrices
    mass_matrix = new double[col * col];
    memset(mass_matrix, 0, sizeof(double) * col * col);
    gradC_matrix = new double[row * col];
    memset(gradC_matrix, 0, sizeof(double) * row * col);
    
    gradCDot_matrix = new double[row * col];
    memset(gradCDot_matrix, 0, sizeof(double) * row * col);
    
    trans_gradC_matrix = new double[col * row];
    memset(trans_gradC_matrix, 0, sizeof(double) * col * row);
    
    int new_dimension = col+row;
    left_matrix = new double[new_dimension * new_dimension];
    memset(left_matrix, 0, sizeof(double)*new_dimension*new_dimension);
    
    ext_force_matrix = new double[col];
    
    qDot_matrix = new double[col];
    memset(qDot_matrix, 0, sizeof(double) * col);
    
    gradCq_DotDot_matrix1 = new double[row];
    memset(gradCq_DotDot_matrix1, 0, sizeof(double) * row);
    
    gradCq_DotDot_matrix2 = new double[row];
    memset(gradCq_DotDot_matrix2, 0, sizeof(double) * row);
    
    gradCq_DotDot_matrix3 = new double[row];
    memset(gradCq_DotDot_matrix3, 0, sizeof(double) * row);
    
    right_matrix = new double[new_dimension];
    memset(right_matrix, 0, sizeof(double) * new_dimension);
    
    x = gsl_vector_alloc(new_dimension);
    gsl_vector_set_zero(x);
};

ParticleSystem::~ParticleSystem()
{
    delete [] mass_matrix;
    delete [] gradC_matrix;
    delete [] gradCDot_matrix;
    delete [] trans_gradC_matrix;
    delete [] left_matrix;
    delete [] ext_force_matrix;
    delete [] qDot_matrix;
    delete [] gradCq_DotDot_matrix1;
    delete [] gradCq_DotDot_matrix2;
    delete [] gradCq_DotDot_matrix3;
    delete [] right_matrix;
    gsl_vector_free(x);
};

double ParticleSystem::linkConstraint(int idx)
{
    // from the third particle in the list(idx = 2), particle_lst[0].pnt_position = (0, 0, 0)T
    double x = particle_lst[idx - 1].pnt_position.x - particle_lst[idx].pnt_position.x;
    double y = particle_lst[idx - 1].pnt_position.y - particle_lst[idx].pnt_position.y;
    double c  = x*x + y*y - len_link*len_link;
    return c;
};

double ParticleSystem::ringConstraint(int idx)
{
    double c = particle_lst[idx].pnt_position.x * particle_lst[idx].pnt_position.x + particle_lst[idx].pnt_position.y * particle_lst[idx].pnt_position.y + particle_lst[idx].pnt_position.y;
    return c;
};

void ParticleSystem::linkGradient(int idx)
{
    link_c_q[0] = 2 * (particle_lst[idx-1].pnt_position.x - particle_lst[idx].pnt_position.x);
    link_c_q[1] = 2 * (particle_lst[idx-1].pnt_position.y - particle_lst[idx].pnt_position.y);
    link_c_q[2] = 2 * (particle_lst[idx].pnt_position.x - particle_lst[idx-1].pnt_position.x);
    link_c_q[3] = 2 * (particle_lst[idx].pnt_position.y - particle_lst[idx-1].pnt_position.y);
};

void ParticleSystem::ringGradient(int idx)
{
    ring_c_q[0] = 2 * particle_lst[idx].pnt_position.x;
    ring_c_q[1] = 2 * particle_lst[idx].pnt_position.y + 1;
};

void ParticleSystem::computeLinkCDot(int idx)
{
    link_cDot_q[0] = 2 * (particle_lst[idx-1].pnt_velocity.x - particle_lst[idx].pnt_velocity.x);
    link_cDot_q[1] = 2 * (particle_lst[idx-1].pnt_velocity.y - particle_lst[idx].pnt_velocity.y);
    link_cDot_q[2] = 2 * (particle_lst[idx].pnt_velocity.x - particle_lst[idx-1].pnt_velocity.x);
    link_cDot_q[3] = 2 * (particle_lst[idx].pnt_velocity.y - particle_lst[idx-1].pnt_velocity.y);
};

void ParticleSystem::computeRingCDot(int idx)
{
    ring_cDot_q[0] = 2 * particle_lst[idx].pnt_velocity.x;
    ring_cDot_q[1] = 2 * particle_lst[idx].pnt_velocity.y;
};

void ParticleSystem::massMatrix()
{
    double m = particle_lst[1].mass_pnt;
    for (int i = 0, j = 0; i < col; i++,j++)
    {
        mass_matrix[i*col+j] = m;
    }
};

void ParticleSystem::gradCMatrix()
{
    int idx = 0;
    int idx_n = no_particle - 1; // the last in the particle list
    int dim = row*col; // gradientC matrix is a row * col matrix
    // the second
    gradC_matrix[0] = 2 * particle_lst[1].pnt_position.x;
    gradC_matrix[1] = 2 * particle_lst[1].pnt_position.y;
    // from the third to the last
    for (int i=2; i<no_particle; i++)
    {
        idx = (i-1) * col + 2 * (i-2);
        linkGradient(i);
        gradC_matrix[idx] = link_c_q[0];
        gradC_matrix[idx+1] = link_c_q[1];
        gradC_matrix[idx+2] = link_c_q[2];
        gradC_matrix[idx+3] = link_c_q[3];
    }
    // ring constraint to the last
    ringGradient(idx_n);
    gradC_matrix[dim-2] = ring_c_q[0];
    gradC_matrix[dim-1] = ring_c_q[1];
};

void ParticleSystem::gradCDotMatrix()
{
    int idx = 0;
    int idx_n = no_particle-1; // the last particle in the particle list
    int dim = row * col; // gradientC matrix is a row * col matrix
    // the first particle
    gradCDot_matrix[0] = 2 * particle_lst[1].pnt_velocity.x;
    gradCDot_matrix[1] = 2 * particle_lst[1].pnt_velocity.y;
    // from the second to the last
    for (int i = 2; i < no_particle; i++)
    {
        idx = (i-1) * col + 2 * (i-2);
        
        computeLinkCDot(i);
        gradCDot_matrix[idx] = link_cDot_q[0];
        gradCDot_matrix[idx + 1] = link_cDot_q[1];
        gradCDot_matrix[idx + 2] = link_cDot_q[2];
        gradCDot_matrix[idx + 3] = link_cDot_q[3];
    }
    // ring constraint to the last
    computeRingCDot(idx_n);
    gradCDot_matrix[dim - 2] = ring_cDot_q[0];
    gradCDot_matrix[dim - 1] = ring_cDot_q[1];
};

void ParticleSystem::transGradCMatrix()
{
    for (int i=0; i<row; i++)
    {
        for (int j=0; j<col; j++)
        {
            trans_gradC_matrix[j*row+i] = gradC_matrix[i*col+j];
        }
    }
};

void ParticleSystem::combineToLeft()
{
    int new_dimension = row + col;
    
    // insert mass matrix to left matrix
    massMatrix();
    int i, j, k;
    for (i=0; i<col; i++)
    {
        for (int j=0; j<col; j++)
        {
            left_matrix[i*new_dimension+j] = mass_matrix[i*col+j];
        }
    }
    
    // insert gradientC matrix
    gradCMatrix();
    
    for (i=col, k=0; i<new_dimension&& k<row; i++, k++)
    {
        for (int j=0; j<col; j++)
        {
            left_matrix[i*new_dimension+j] = gradC_matrix[k*col+j];
        }
    }
    
    // insert transposed gradientC matrix
    transGradCMatrix();
    for (int i=0; i<col; i++)
    {
        for (j=col, k=0; j<new_dimension&& k<row; j++, k++)
        {
            left_matrix[i*new_dimension+j] = trans_gradC_matrix[i*row+k];
        }
    }
};

void ParticleSystem::extForceMatrix()
{
    memset(ext_force_matrix, 0, sizeof(double) * col);
    double m = particle_lst[1].mass_pnt;
    
    for (int i=1; i<col; i+=2)
    {
        ext_force_matrix[i] = ext_force_matrix[i] + gravity.y*m;
    }
    
    // UI control
    for (int i=0; i<col; i+=2)
    {
        ext_force_matrix[i] = ext_force_matrix[i]-left_ext_force*m;
        ext_force_matrix[i] = ext_force_matrix[i]+right_ext_force*m;
        ext_force_matrix[i+1] = ext_force_matrix[i+1]+up_ext_force*m;
        ext_force_matrix[i+1] = ext_force_matrix[i+1]-down_ext_force*m;
    }
};

void ParticleSystem::matrixMultiple(double *mat_top, double *mat_new, double *mat_out)
{
    // allocate memory
    gsl_matrix *gsl_top = gsl_matrix_alloc(row, col);
    gsl_vector *gsl_new = gsl_vector_alloc(col);
    gsl_vector *gsl_out = gsl_vector_alloc(row);
    
    // set top
    for (int i=0; i<row; i++)
    {
        for (int j=0; j<col; j++)
        {
            gsl_matrix_set(gsl_top, i, j, mat_top[i*col+j]);
        }
    }
    
    // set new
    for (int i=0; i<col; i++)
    {
        gsl_vector_set(gsl_new, i, mat_new[i]);
    }
    
    // multiplication
    gsl_blas_dgemv(CblasNoTrans, 1, gsl_top, gsl_new, 0.0, gsl_out);
    //	gsl_blas_dgemv(CblasNoTrans, 1, gsl_top, gsl_new, 0.0, gsl_out);
    for (int i=0; i<row; i++)
    {
        mat_out[i] = gsl_vector_get(gsl_out, i);
    }
    gsl_matrix_free(gsl_top);
    gsl_vector_free(gsl_new);
    gsl_vector_free(gsl_out);
};

void ParticleSystem::computeRight()
{
    int new_dimension = row + col;
    int i, j;
    int idx_n = no_particle - 1;
    // the last in the particle list
    
    extForceMatrix();
    
    // insert external force
    for(int i=0; i<col; i++)
    {
        right_matrix[i] = ext_force_matrix[i];
    }
    // matrix q'
    for(i=0, j=1; i<col&& j<no_particle; i+=2, j++)
    {
        qDot_matrix[i] = particle_lst[j].pnt_velocity.x;
        qDot_matrix[i+1] = particle_lst[j].pnt_velocity.y;
    }
    // matrix dC' / dq
    gradCDotMatrix();
    //-(dC' / dq)q'
    matrixMultiple(gradCDot_matrix, qDot_matrix, gradCq_DotDot_matrix1);
    for(i=0; i<row; i++)
    {
        gradCq_DotDot_matrix1[i] = (-1.0) * gradCq_DotDot_matrix1[i];
    }
    // compute (dC / dq)q'
    matrixMultiple(gradC_matrix, qDot_matrix, gradCq_DotDot_matrix2);
    // compute C
    gradCq_DotDot_matrix3[0] = particle_lst[1].pnt_position.x * particle_lst[1].pnt_position.x + particle_lst[1].pnt_position.y * particle_lst[1].pnt_position.y - len_link * len_link;
    for(i=2; i<no_particle; i++)
    {
        gradCq_DotDot_matrix3[i-1] = linkConstraint(i);
    }
    gradCq_DotDot_matrix3[row-1] = ringConstraint(idx_n);
    // insert -(dC' / dq)q' - alpha * (dC / dq)q' - beta * C
    for(i=0; i<row; i++)
    {
        right_matrix[i+col] = gradCq_DotDot_matrix1[i] - alpha * gradCq_DotDot_matrix2[i] - beta * gradCq_DotDot_matrix3[i];
    }
};

void ParticleSystem::computeAcceleration()
{
    int new_dimension = row + col;
    combineToLeft();
    computeRight();
    // allocate memory
    gsl_matrix *gsl_v = gsl_matrix_alloc(new_dimension, new_dimension);
    gsl_vector *gsl_s = gsl_vector_alloc(new_dimension);
    gsl_vector *gsl_work = gsl_vector_alloc(new_dimension);
    gsl_matrix_view m = gsl_matrix_view_array(left_matrix, new_dimension, new_dimension);
    gsl_vector_view b = gsl_vector_view_array(right_matrix, new_dimension);
    gsl_linalg_SV_decomp(&m.matrix, gsl_v, gsl_s, gsl_work);
    
    for (int i=0; i<new_dimension; i++)
    {
        if (fabs(gsl_vector_get(gsl_s, i) / gsl_vector_get(gsl_s, 0)) < 1E-6)
        {
            gsl_vector_set(gsl_s, i, 0.0);
        }
    }
    gsl_linalg_SV_solve(&m.matrix, gsl_v, gsl_s, &b.vector, x);
    
    for (int i=0, j=1; i<col&& j<no_particle; i+=2, j++)
    {
        particle_lst[j].pnt_acceleration.x = gsl_vector_get(x, i);
        particle_lst[j].pnt_acceleration.y = gsl_vector_get(x, i + 1);
    }
    
    gsl_matrix_free(gsl_v);
    gsl_vector_free(gsl_s);
    gsl_vector_free(gsl_work);
    
    // reset
    memset(left_matrix, 0, sizeof(double)*new_dimension*new_dimension);
    memset(right_matrix, 0, sizeof(double)*new_dimension);
    gsl_vector_set_zero(x);
};

/*void ParticleSystem::PrintError()
 {
	double error = 0.0;
	for (int i = 0; i < no_particle; ++i)
	{
 // constraint may not always be zero gradCq_DotDot_matrix1
 error += gradCq_DotDot_matrix3[i];
 }
	error_value << error << endl;
	frame++;
	if (frame >= 300)
	{
 error_value.close();
	}
 }*/
void ParticleSystem::glDrawCylin(float c_radius, int c_subdivision)
{
    float vx, vy, vz, v, ax;
    
    GLUquadricObj *quadratic = gluNewQuadric();
    gluQuadricNormals(quadratic, GLU_SMOOTH);
    
    for(int i=0, j=1; j<no_particle; i++, j++)
    {
        vx = particle_lst[j].pnt_position.x - particle_lst[i].pnt_position.x;
        vy = particle_lst[j].pnt_position.y - particle_lst[i].pnt_position.y;
        vz = particle_lst[j].pnt_position.z - particle_lst[i].pnt_position.z;
        v = sqrt(vx*vx + vy*vy + vz*vz);
        
        if (fabs(vz)<1.0e-3)
        {
            ax = 57.2957 * acos(vx/v);
            if (vy <= 0.0)
                ax*=-1;
        }else
        {
            ax = 57.2957 * acos(vz/v);
            if (vz<=0.0)
                ax*=-1;
        }
        
        float rx = -vy * vz;
        float ry = vx * vz;
        glPushMatrix();
        glTranslated(particle_lst[i].pnt_position.x, particle_lst[i].pnt_position.y, particle_lst[i].pnt_position.z);
        if (fabs(vz) < 1.0e-3)
        {
            glRotated(90.0, 0, 1, 0.0);
            glRotated(ax, -1.0, 0.0, 0.0);
        }
        else
        {
            glRotated(ax, rx, ry, 0.0);
        }
        gluQuadricOrientation(quadratic, GLU_OUTSIDE);
        gluCylinder(quadratic, c_radius, c_radius, v, c_subdivision, 1);
        gluQuadricOrientation(quadratic, GLU_INSIDE);
        gluDisk(quadratic, 0.0, c_radius, c_subdivision, 1);
        glTranslated(0, 0, v);
        gluQuadricOrientation(quadratic, GLU_OUTSIDE);
        gluDisk(quadratic, 0.0, c_radius, c_subdivision, 1);
        glPopMatrix();
    }
    gluDeleteQuadric(quadratic);
};

void ParticleSystem::glDrawSphere(float s_radius, int s_subdivisions)
{
    for (int i = 0; i < no_particle; i++)
    {
        GLUquadricObj *quadratic = gluNewQuadric();
        gluQuadricNormals(quadratic, GLU_SMOOTH);
        glPushMatrix();
        glTranslated(particle_lst[i].pnt_position.x, particle_lst[i].pnt_position.y, particle_lst[i].pnt_position.z);
        gluSphere(quadratic, s_radius, s_subdivisions, s_subdivisions);
        glPopMatrix();
        gluDeleteQuadric(quadratic);
    }
    
}

void ParticleSystem::glDrawRing()
{
    point normal;
    double major_step = 2.0f*pi / no_major;
    double minor_step = 2.0f*pi / no_minor;
    for (int i=0; i<no_major; i++)
    {
        double a0 = i * major_step;
        double a1 = a0 + major_step;
        float x0 = (float)cos(a0);
        float y0 = (float)sin(a0);
        float x1 = (float)cos(a1);
        float y1 = (float)sin(a1);
        glBegin(GL_TRIANGLE_STRIP);
        for (int j=0; j<=no_minor; j++)
        {
            double b = j * minor_step;
            float c = (float) cos(b);
            float r = in_radius * c + out_radius;
            float z = in_radius * (float) sin(b);
            
            glTexCoord2f((float)i / (float)(no_major), (float)(j) / (float)(no_minor));
            normal.x = x0 * c;
            normal.y = y0 * c;
            normal.z = z / in_radius;
            double len = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
            normal.x = normal.x / len;
            normal.y = normal.y / len;
            normal.z = normal.z / len;
            
            glNormal3f(normal.x, normal.y, normal.z);
            glVertex3f(x0 * r + ring_pos.x, y0 * r + ring_pos.y, z + ring_pos.z);
            
            glTexCoord2f((float)(i + 1) / (float)(no_major), (float)(j) / (float)(no_minor));
            normal.x = x1 * c;
            normal.y = y1 * c;
            normal.z = z / in_radius;
            
            len = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
            normal.x = normal.x / len;
            normal.y = normal.y / len;
            normal.z = normal.z / len;
            
            glNormal3f(normal.x, normal.y, normal.z);
            glVertex3f(x1 * r + ring_pos.x, y1 * r + ring_pos.y, z + ring_pos.z);
        }
        glEnd();
    }
}

void ParticleSystem::glRender(float s_radius, int s_subdivision, float c_radius, int c_subdivision)
{
    glDrawCylin(c_radius, c_subdivision);
    glDrawSphere(s_radius,s_subdivision);
    glDrawRing();
};



