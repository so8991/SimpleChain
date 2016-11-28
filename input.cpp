/*
 
 USC/Viterbi/Computer Science
 "Jello Cube" Assignment 1 starter code
 
 Your name:
 Sunoh Yoo
 
 */

#include "jello.h"
#include "input.h"
#include "particlesystem.h"
#include "pic.h"

extern ParticleSystem *particle_system;

/* Write a screenshot, in the PPM format, to the specified filename, in PPM format */
void saveScreenshot(int windowWidth, int windowHeight, char *filename)
{
    if (filename == NULL)
        return;
    
    // Allocate a picture buffer
    Pic * in = pic_alloc(windowWidth, windowHeight, 3, NULL);
    
    printf("File to save to: %s\n", filename);
    
    for (int i=windowHeight-1; i>=0; i--)
    {
        glReadPixels(0, windowHeight-i-1, windowWidth, 1, GL_RGB, GL_UNSIGNED_BYTE,
                     &in->pix[i*in->nx*in->bpp]);
    }
    
    if (ppm_write(filename, in))
        printf("File saved Successfully\n");
    else
        printf("Error in Saving\n");
    
    pic_free(in);
}

vector<int> getMouseLocation(int x, int y)
{
    vector<int> location;
    location.push_back(x);
    location.push_back(y);
    return location;
}
/* converts mouse drags into information about rotation/translation/scaling */
void mouseMotionDrag(int x, int y)
{
    printf("x:%d, y:%d\n", x, y);
    int vMouseDelta[2] = {x-g_vMousePos[0], y-g_vMousePos[1]};
    
    if (g_iRightMouseButton) // handle camera rotations
    {
        Phi += vMouseDelta[0] * 0.01;
        Theta += vMouseDelta[1] * 0.01;
        
        if (Phi>2*pi)
            Phi -= 2*pi;
        
        if (Phi<0)
            Phi += 2*pi;
        
        if (Theta>pi / 2 - 0.01) // dont let the point enter the north pole
            Theta = pi / 2 - 0.01;
        
        if (Theta<- pi / 2 + 0.01)
            Theta = -pi / 2 + 0.01;
        
        g_vMousePos[0] = x;
        g_vMousePos[1] = y;
    }
}

void mouseMotion (int x, int y)
{
    g_vMousePos[0] = x;
    g_vMousePos[1] = y;
}

void mouseButton(int button, int state, int x, int y)
{
    switch (button)
    {
        case GLUT_LEFT_BUTTON:
            g_iLeftMouseButton = (state==GLUT_DOWN);
            break;
        case GLUT_MIDDLE_BUTTON:
            g_iMiddleMouseButton = (state==GLUT_DOWN);
            break;
        case GLUT_RIGHT_BUTTON:
            g_iRightMouseButton = (state==GLUT_DOWN);
            break;
    }
    
    g_vMousePos[0] = x;
    g_vMousePos[1] = y;
}


// gets called whenever a key is pressed
void keyboardFunc (unsigned char key, int x, int y)
{
    switch (key)
    {
        case 27:
            exit(0);
            break;
            
        case 'e':
             Theta = pi / 6;
             Phi = pi / 6;
             viewingMode = 0;
             break;
            
            //case 'p':
            //   = 1 - pause;
            //   break;
            
        case 'z':
            R -= 0.2;
            if (R < 0.2)
                R = 0.2;
            break;
            
        case 'x':
            R += 0.2;
            break;
            
        case ' ':
            saveScreenToFile = 1 - saveScreenToFile;
            break;
            
            // reset external force
        case 'r':
            particle_system->up_ext_force = 0.0;
            particle_system->down_ext_force = 0.0;
            particle_system->left_ext_force = 0.0;
            particle_system->right_ext_force = 0.0;
            break;
            
            // switch 'Euler' to 'RK4'
        case 'k':
            particle_system->integrator = 'R';
            break;
            
            // switch 'RK4' to 'Euler'
        case 'u':
            particle_system->integrator = 'E';
            break;
    }
}

void SpecialKeys(int key, int x, int y)
{
    if (key == GLUT_KEY_UP)
    {
        particle_system->up_ext_force += 0.1;
    }
    
    if (key == GLUT_KEY_DOWN)
    {
        particle_system->down_ext_force += 0.1;
    }
    
    if (key == GLUT_KEY_LEFT)
    {
        particle_system->left_ext_force += 0.1;
    }
    
    if (key == GLUT_KEY_RIGHT)
    {
        particle_system->right_ext_force += 0.1;
    }
}

