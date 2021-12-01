#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdio.h>
#include "glew.h"
#include "cuda_gl_interop.h"
#include "freeglut.h"


static int windowID;
static float width, height;
// static int widthINT, heightINT;
__device__ float *arrayCUDA;
float *arrayHOST;
int *parrayCUDA;
float *parrayHOST;
static int arraySize, arrayDimension;

static float dt, visc;

static void Sliders()
{
    // Drawing Sliders Text Fields
    glColor3f(0.0, 0.0, 0.0);
    glRasterPos2f(408., 26.);
    glutBitmapString(GLUT_BITMAP_9_BY_15, "Controls");
    glRasterPos2f(350., 48.);
    glutBitmapString(GLUT_BITMAP_8_BY_13, "Time Step");
    glRasterPos2f(350., 68.);
    glutBitmapString(GLUT_BITMAP_8_BY_13, "Viscosity");
    glRasterPos2f(350., 88.);
    glutBitmapString(GLUT_BITMAP_8_BY_13, "Vorticity");
    glRasterPos2f(0., 0.);

    glBegin(GL_LINES);

    // Time Step
    glVertex2d(430., 38.);
    glVertex2d(510., 38.);
    glVertex2d(430., 49.);
    glVertex2d(510., 49.);
    glVertex2d(430., 38.);
    glVertex2d(430., 49.);
    glVertex2d(510., 38.);
    glVertex2d(510., 49.);

    // Viscosity
    glVertex2d(430., 58.);
    glVertex2d(510., 58.);
    glVertex2d(430., 69.);
    glVertex2d(510., 69.);
    glVertex2d(430., 58.);
    glVertex2d(430., 69.);
    glVertex2d(510., 58.);
    glVertex2d(510., 69.);

    // Vorticity
    glVertex2d(430., 78.);
    glVertex2d(510., 78.);
    glVertex2d(430., 89.);
    glVertex2d(510., 89.);
    glVertex2d(430., 78.);
    glVertex2d(430., 89.);
    glVertex2d(510., 78.);
    glVertex2d(510., 89.);

    // Fill In Sliders
    float sliderStart = 430.;
    float sliderEnd = 510.;
    if (dt > 30.)
        dt = 30.;
    if (visc > 1.)
        visc = 1.;
    float dtSliderEnd = ((dt / 30.) * 80.) + sliderStart;
    float viscSliderEnd = ((visc / 1.) * 80.) + sliderStart;
    for (float i = sliderStart; i <= sliderEnd; i++)
    {
        if (i <= dtSliderEnd)
        {
            glVertex2d(i, 38.);
            glVertex2d(i, 49.);
        }
        if (i <= viscSliderEnd)
        {
            glVertex2d(i, 58.);
            glVertex2d(i, 69.);
        }
    }
    glEnd();
}

__global__ static void color_Array()
{
    // int index = blockIdx.x * blockDim.x + threadIdx.x;
    for (int i = 0; i < 512; i++)
    {
        for (int j = 0; j < 512; j++)
        {
            arrayCUDA[i + j * 512] = 0.5;
        }
    }
}
static void draw_Array()
{
    glBegin(GL_POINTS);

    /////// This is where the problem lies
    // How to copy a __device__ array to the __host__ array.
    cudaMemcpyFromSymbol(arrayHOST, (const char *)"arrayCUDA", arraySize, 0, cudaMemcpyDeviceToHost);
    ///////
    for (int i = 0; i < 512; i++)
    {
        for (int j = 0; j < 512; j++)
        {
            float color = arrayHOST[i + j * 512];
            glColor3f(color, color, color);
            glVertex2i(i, j);
        }
    }
    glEnd();
}

static void pre_display(void)
{
    glViewport(0, 0, 512, 512);
    glLoadIdentity();
    gluOrtho2D(0.0, 512., 512., 0.0);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
}

static void post_display(void)
{
    glutSwapBuffers();
}


static void key_func(unsigned char key, int x, int y)
{
    switch (key)
    {
    case 'w':
    case 'W':
        dt += 0.1;
        printf("dt is now %f\n", dt);
        break;

    case 's':
    case 'S':
        dt -= 0.1;
        printf("dt is now %f\n", dt);
        break;

    case 'e':
    case 'E':
        visc += 0.05;
        printf("visc is now %f\n", visc);
        break;

    case 'd':
    case 'D':
        visc -= 0.05;
        printf("visc is now %f\n", visc);
        break;
    }
}

static void idle_func(void)
{
    glutSetWindow(windowID);
    glutPostRedisplay();
}

static void display_func(void)
{
    pre_display();
    color_Array<<<1,1>>>();
    cudaDeviceSynchronize();
    draw_Array();
    Sliders();
    post_display();
}

static void open_glut_window(void)
{
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(512, 512);
    windowID = glutCreateWindow("Fluid Simulation");
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glutSwapBuffers();
    glClear(GL_COLOR_BUFFER_BIT);
    glutSwapBuffers();
    pre_display();
    glutKeyboardFunc(key_func);
    glutIdleFunc(idle_func);
    glutDisplayFunc(display_func);
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    width = 512.;
    height = 512.;
    dt = 0.1;
    visc = 0.;
    arrayDimension = width * height;
    arraySize = arrayDimension * sizeof(float);
    cudaMalloc(&arrayCUDA, arraySize);
    arrayHOST = (float *)malloc(arraySize);
    std::cout << "arrayCUDA " << &arrayCUDA << "\n";
    std::cout << "arrayCUDAP " << parrayCUDA << "\n";
    std::cout << "arrayCUDAPP " << &parrayCUDA << "\n";
    open_glut_window();
    glutMainLoop();
    // free(arrayCUDA);
    free(arrayHOST);
    exit(0);
}
