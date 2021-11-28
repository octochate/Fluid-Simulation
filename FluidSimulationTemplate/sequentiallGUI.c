#include <stdlib.h>
#include <stdio.h>
#include "freeglut.h"

static int windowID;
static float width, height;
static int widthINT, heightINT;
static float *array;
static int arraySize, arrayDimension;

static float dt, visc;

static void RenderString()
{
    // glColor3f(1.0, 1.0, 1.0);
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
    glVertex2d(430., 38.);
    glVertex2d(510., 38.);
    glVertex2d(430., 49.);
    glVertex2d(510., 49.);
    glVertex2d(430., 38.);
    glVertex2d(430., 49.);
    glVertex2d(510., 38.);
    glVertex2d(510., 49.);
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

    // First
    glVertex2d(430., 58.);
    glVertex2d(510., 58.);
    glVertex2d(430., 69.);
    glVertex2d(510., 69.);
    glVertex2d(430., 58.);
    glVertex2d(430., 69.);
    glVertex2d(510., 58.);
    glVertex2d(510., 69.);
    // First
    glVertex2d(430., 78.);
    glVertex2d(510., 78.);
    glVertex2d(430., 89.);
    glVertex2d(510., 89.);
    glVertex2d(430., 78.);
    glVertex2d(430., 89.);
    glVertex2d(510., 78.);
    glVertex2d(510., 89.);

    glEnd();
}

static void color_Array()
{
    // glBegin(GL_QUADS);
    // glBegin(GL_POINTS);
    for (int i = 0; i < 512; i++)
    {
        for (int j = 0; j < 512; j++)
        {
            // float color = (float) (i * j) / (float) arrayDimension;
            // glColor3f(1.0, 0.0, 0.0);
            // glColor3f(color, color, color);
            array[i + j * 512] = 255;
            // float color1 = rand() % 255;
            // float color2 = rand() % 255;
            // float color3 = rand() % 255;
            // glColor3f(color1, color2, color3);
            // glColor3i(color, color, color);
            // glVertex2i(i, j);

            // glVertex2d(0., 0.);
            // glVertex2d(256., 256.);
        }
    }
    // glEnd();
}
static void draw_Array()
{
    // glBegin(GL_QUADS);
    glBegin(GL_POINTS);

    for (int i = 0; i < 512; i++)
    {
        for (int j = 0; j < 512; j++)
        {
            float color = array[i + j * 512];
            glColor3f(color, color, color);
            glVertex2i(i, j);
        }
    }
    glEnd();
}

static void pre_display(void)
{
    glViewport(0, 0, 512, 512);
    // glMatrixMode(GL_PROJECTION);
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
    color_Array();
    glutPostRedisplay();
}

static void display_func(void)
{
    pre_display();
    draw_Array();
    RenderString();
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
    // cudaMallocManaged(&array, width * height * sizeof(double));
    array = (float *)malloc(arraySize);
    open_glut_window();
    glutMainLoop();
    exit(0);
}