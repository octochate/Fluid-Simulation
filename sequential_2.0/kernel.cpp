<<<<<<< HEAD
﻿#include <stdlib.h>
=======
#include <stdlib.h>
>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "freeglut.h"
#include "utilities.h"


<<<<<<< HEAD


=======
>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62
#define NUM_ROW 64
#define NUM_COL NUM_ROW
#define GRID_SIZE NUM_ROW * NUM_COL
#define BOUNDARY_GAIN 0.75

volatile int q = 0;
volatile int q_ink = 0;
<<<<<<< HEAD
=======

>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62
volatile int dvel;

static int win_id;
static int win_x, win_y;
static float dt, diff, visc;

<<<<<<< HEAD
=======
static int diffusion_iter = 50, pressure_iter = 50;

>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62
static float vort, fluidSize, fluidAmount, forceSize, forceAmount;

// Variable bounds.
float numIntervals = 100;
float dtMax = 0.5;
float dtMin = 0.000001;
float dt_del = (dtMax - dtMin) / numIntervals;
<<<<<<< HEAD
float viscMax = 0.9;
float viscMin = 0.00001;
=======
float viscMax = 1;
float viscMin = 0.000000001;
>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62
float visc_del = (viscMax - viscMin) / numIntervals;
float vortMax = 1.;
float vortMin = 0.001;
float vort_del = (vortMax - vortMin) / numIntervals;
float fluidSizeMax = 5.;
float fluidSizeMin = 0.01;
float fluidSize_del = (fluidSizeMax - fluidSizeMin) / numIntervals;
float fluidAmountMax = 1.;
float fluidAmountMin = 0.01;
float fluidAmount_del = (fluidAmountMax - fluidAmountMin) / numIntervals;
float forceSizeMax = 5.;
float forceSizeMin = 0.01;
float forceSize_del = (forceSizeMax - forceSizeMin) / numIntervals;
float forceAmountMax = 1.;
float forceAmountMin = 0.01;
float forceAmount_del = (forceAmountMax - forceAmountMin) / numIntervals;

static void get_from_UI();
static void display_func(void);

// update functions
<<<<<<< HEAD
=======
void advect_density(f3 u[NUM_ROW][NUM_COL], f3 u1[NUM_ROW][NUM_COL], f3 X[NUM_ROW][NUM_COL], float timestep, float gain);
>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62
void update_grid(f3 u[NUM_ROW][NUM_COL], f3 u_temp[NUM_ROW][NUM_COL], f3 p[NUM_ROW][NUM_COL], f3 divergence_u[NUM_ROW][NUM_COL], float timestep, float epsilon);
void advect(f3 u[NUM_ROW][NUM_COL], f3 u1[NUM_ROW][NUM_COL], f3 X[NUM_ROW][NUM_COL], float timestep, float gain);
void diffuse(f3 u[NUM_ROW][NUM_COL], f3 u1[NUM_ROW][NUM_COL], float diffusion_rate, float timestep);
void addForces(f3 u[NUM_ROW][NUM_COL], float timestep);
void computePressure(f3 p[NUM_ROW][NUM_COL], f3 u[NUM_ROW][NUM_COL], f3 divergence_u[NUM_ROW][NUM_COL], float timestep);
void subtractPressureGradient(f3 u[NUM_ROW][NUM_COL], f3 u1[NUM_ROW][NUM_COL], f3 p[NUM_ROW][NUM_COL]);
void divergence(f3 div_f[NUM_ROW][NUM_COL], f3 f[NUM_ROW][NUM_COL]);

// GUI event
void addExternalForce(f3 force, f3 position, float impulseRadius);
void addInk(f3 force, f3 position, float impulseRadius);

// assign pointers for u, u1, u2 2D grids
<<<<<<< HEAD
f3 U[2][NUM_ROW][NUM_COL], P[NUM_ROW][NUM_COL], P_guess[NUM_ROW][NUM_COL], Divergence_u[NUM_ROW][NUM_COL], ink[2][NUM_ROW][NUM_COL];
=======
f3 U[2][NUM_ROW][NUM_COL], U_b[NUM_ROW][NUM_COL], P[NUM_ROW][NUM_COL], P_guess[NUM_ROW][NUM_COL], Divergence_u[NUM_ROW][NUM_COL], ink[2][NUM_ROW][NUM_COL];
>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62

void clearGrids() {
    for (int x = 0; x < NUM_COL; x++) {
        for (int y = 0; y < NUM_ROW; y++) {
            P[y][x] = Divergence_u[y][x] = U[0][y][x] = U[1][y][x] = ink[0][y][x] = ink[1][y][x] = f3(0,0);
        }
    }
}

void RenderString(const char* string, float x, float y)
{
	glColor3f(1.0, 1.0, 1.0);
	glRasterPos2f(x, y);
}

static void Sliders()
{
    const unsigned char a[50] = "Controls";
    const unsigned char b[50] = "Time Step";
    const unsigned char c[50] = "Viscosity";
    const unsigned char d[50] = "Vorticity";
    const unsigned char e[50] = "Fluid Size";
    const unsigned char f[50] = "Fluid Amount";
    const unsigned char g[50] = "Force Size";
    const unsigned char h[50] = "Force Amount";

    const unsigned char* aPtr = a;
    const unsigned char* bPtr = b;
    const unsigned char* cPtr = c;
    const unsigned char* dPtr = d;
    const unsigned char* ePtr = e;
    const unsigned char* fPtr = f;
    const unsigned char* gPtr = g;
    const unsigned char* hPtr = h;

    // Drawing Sliders Text Fields
    glColor3f(1.0, 1.0, 1.0);
    glRasterPos2f(0.78125, 0.94921875);
    glutBitmapString(GLUT_BITMAP_9_BY_15, aPtr);
    glRasterPos2f(0.64453125, 0.90625000);
    glutBitmapString(GLUT_BITMAP_8_BY_13, bPtr);
    glRasterPos2f(0.64453125, 0.86718750);
    glutBitmapString(GLUT_BITMAP_8_BY_13, cPtr);
    glRasterPos2f(0.64453125, 0.82812500);
    glutBitmapString(GLUT_BITMAP_8_BY_13, dPtr);
    glRasterPos2f(0.64453125, 0.7890625);
    glutBitmapString(GLUT_BITMAP_8_BY_13, ePtr);
    glRasterPos2f(0.64453125, 0.7500000);
    glutBitmapString(GLUT_BITMAP_8_BY_13, fPtr);
    glRasterPos2f(0.64453125, 0.7109375);
    glutBitmapString(GLUT_BITMAP_8_BY_13, gPtr);
    glRasterPos2f(0.64453125, 0.6718750);
    glutBitmapString(GLUT_BITMAP_8_BY_13, hPtr);

    glRasterPos2f(0., 0.);

    glBegin(GL_LINES);
    glColor3f(1.0, 1.0, 1.0);

    // Draw slider boxes.
    for (int i = 0; i < 7; i++)
    {
        // Compute heights.
        float heightTop = 1. - (38. + (float)i * 20.) / 512.;
        float heightBottom = 1. - (49. + (float)i * 20.) / 512.;
        glVertex2d(0.83984375, heightTop);
        glVertex2d(0.99609375, heightTop);
        glVertex2d(0.83984375, heightBottom);
        glVertex2d(0.99609375, heightBottom);
        glVertex2d(0.83984375, heightTop);
        glVertex2d(0.83984375, heightBottom);
        glVertex2d(0.99609375, heightTop);
        glVertex2d(0.99609375, heightBottom);
    }


    // Fill In Sliders
    float sliderStart = 0.83984375;
    float sliderEnd = 0.99609375;
    // Variable bounds.

    // Stay within bounds.
    if (dt > dtMax)
        dt = dtMax;
    if (dt < dtMin)
        dt = dtMin;
    if (visc > viscMax)
        visc = viscMax;
    if (visc < viscMin)
        visc = viscMin;

    // Compute dynamic slider fill.
    float dtSliderEnd = ((dt / dtMax) * 0.15625) + sliderStart;
    float viscSliderEnd = ((visc / viscMax) * 0.15625) + sliderStart;
    float vortSliderEnd = ((vort / vortMax) * 0.15625) + sliderStart;
    float fluidSizeSliderEnd = ((fluidSize / fluidSizeMax) * 0.15625) + sliderStart;
    float fluidAmountSliderEnd = ((fluidAmount / fluidAmountMax) * 0.15625) + sliderStart;
    float forceSizeSliderEnd = ((forceSize / forceSizeMax) * 0.15625) + sliderStart;
    float forceAmountSliderEnd = ((forceAmount / forceAmountMax) * 0.15625) + sliderStart;

    for (float i = sliderStart; i <= sliderEnd; i += 0.001)
    {
        float heightTop = 0.0;
        float heightBottom = 0.0;
        if (i <= dtSliderEnd)
        {
            heightTop = 1. - (38. + 0. * 20.) / 512.;
            heightBottom = 1. - (49. + 0. * 20.) / 512.;
            glVertex2d(i, heightTop);
            glVertex2d(i, heightBottom);
        }
        if (i <= viscSliderEnd)
        {
            heightTop = 1. - (38. + 1. * 20.) / 512.;
            heightBottom = 1. - (49. + 1. * 20.) / 512.;
            glVertex2d(i, heightTop);
            glVertex2d(i, heightBottom);
        }
        if (i <= vortSliderEnd)
        {
            heightTop = 1. - (38. + 2. * 20.) / 512.;
            heightBottom = 1. - (49. + 2. * 20.) / 512.;
            glVertex2d(i, heightTop);
            glVertex2d(i, heightBottom);
        }
        if (i <= fluidSizeSliderEnd)
        {
            heightTop = 1. - (38. + 3. * 20.) / 512.;
            heightBottom = 1. - (49. + 3. * 20.) / 512.;
            glVertex2d(i, heightTop);
            glVertex2d(i, heightBottom);
        }
        if (i <= fluidAmountSliderEnd)
        {
            heightTop = 1. - (38. + 4. * 20.) / 512.;
            heightBottom = 1. - (49. + 4. * 20.) / 512.;
            glVertex2d(i, heightTop);
            glVertex2d(i, heightBottom);
        }
        if (i <= forceSizeSliderEnd)
        {
            heightTop = 1. - (38. + 5. * 20.) / 512.;
            heightBottom = 1. - (49. + 5. * 20.) / 512.;
            glVertex2d(i, heightTop);
            glVertex2d(i, heightBottom);
        }
        if (i <= forceAmountSliderEnd)
        {
            heightTop = 1. - (38. + 6. * 20.) / 512.;
            heightBottom = 1. - (49. + 6. * 20.) / 512.;
            glVertex2d(i, heightTop);
            glVertex2d(i, heightBottom);
        }
    }
    glEnd();
}

/*
  ----------------------------------------------------------------------
   OpenGL specific drawing routines
  ----------------------------------------------------------------------
*/

static void pre_display(void)
{
	glViewport(0, 0, win_x, win_y);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, 1.0, 0.0, 1.0);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
}

static void post_display(void)
{
	glutSwapBuffers();
}

static void draw_velocity(void)
{
	int i, j;
	float x, y, h;

	h = 1.0f / NUM_COL;

	glColor3f(1.0f, 1.0f, 1.0f);
	glLineWidth(1.0f);

	glBegin(GL_LINES);

	for (i = 1; i < NUM_COL; i++) {
		x = (i - 0.5f) * h;
		for (j = 1; j < NUM_ROW; j++) {
			y = (j - 0.5f) * h;

			glVertex2f(x, y);
			glVertex2f(x + U[q % 1][j][i].x, y + U[q % 1][j][i].y);
            //glVertex2f(x + ink[(q_ink) % 1][j][i].x, y + ink[(q_ink) % 1][j][i].y);
		}
	}

	glEnd();
}

int square = 0;

static void draw_density(void)
{
    int i, j;
    float x, y, h;
    f3 d00, d01, d10, d11;

    h = 1.0f / NUM_COL;

    glBegin(GL_QUADS);

    for (i = 0; i < NUM_COL - 1; i++) {
        x = (i - 0.5f) * h;
        for (j = 0; j < NUM_ROW - 1; j++) {
            y = (j - 0.5f) * h;


            if (square == 1 && i >= 34 && i < 62 && j >= 34 && j < 62) {
                ink[(q_ink) % 1][j][i] = f3(0.0, 0);
<<<<<<< HEAD
                d00 = f3(1.0, 10);
                d01 = f3(1.0, 10);
                d10 = f3(1.0, 10);
                d11 = f3(1.0, 10);
                glColor3f(0.0, d00.y, 0.0); glVertex2f(x, y);
                glColor3f(0.0, d10.y, 0.0); glVertex2f(x + h, y);
                glColor3f(0.0, d11.y, 0.0); glVertex2f(x + h, y + h);
                glColor3f(0.0, d01.y, 0.0); glVertex2f(x, y + h);
=======
                d00 = f3(0, 00);
                d01 = f3(0, 00);
                d10 = f3(0, 00);
                d11 = f3(0, 00);
                glColor3f(0.0, d00.x, 0.0); glVertex2f(x, y);
                glColor3f(0.0, d10.x, 0.0); glVertex2f(x + h, y);
                glColor3f(0.0, d11.x, 0.0); glVertex2f(x + h, y + h);
                glColor3f(0.0, d01.x, 0.0); glVertex2f(x, y + h);
>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62
            }
            else {
                d00 = ink[(q_ink) % 1][j][i];
                d01 = ink[(q_ink) % 1][j + 1][i];
                d10 = ink[(q_ink) % 1][j][i + 1];
                d11 = ink[(q_ink) % 1][j + 1][i + 1];

<<<<<<< HEAD
=======
                /*glColor3f(d00.x, d00.y, d00.z); glVertex2f(x, y);
                glColor3f(d10.x, d10.y, d10.z); glVertex2f(x + h, y);
                glColor3f(d11.x, d11.y, d11.z); glVertex2f(x + h, y + h);
                glColor3f(d01.x, d01.y, d01.z); glVertex2f(x, y + h);*/

>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62
                glColor3f(d00.x, d00.y, d00.z); glVertex2f(x, y);
                glColor3f(d10.x, d10.y, d10.z); glVertex2f(x + h, y);
                glColor3f(d11.x, d11.y, d11.z); glVertex2f(x + h, y + h);
                glColor3f(d01.x, d01.y, d01.z); glVertex2f(x, y + h);
            }
        }
    }

    glEnd();
}

/*
  ----------------------------------------------------------------------
   GLUT callback routines
  ----------------------------------------------------------------------
*/


static void reshape_func(int width, int height)
{
	glutSetWindow(win_id);
	glutReshapeWindow(width, height);

	win_x = width;
	win_y = height;
}

static void idle_func(void)
{
    get_from_UI();
    update_grid(U[q % 1], U[(q + 1) % 1], P, Divergence_u, dt, visc);
    q = (q + 1) % 100;


	glutSetWindow(win_id);
	glutPostRedisplay();
}

static void display_func(void)
{
	pre_display();

    if (dvel)
	   draw_velocity();
    else
        draw_density();
    Sliders();
	post_display();
}

/*
  ----------------------------------------------------------------------
   relates mouse movements to forces sources
  ----------------------------------------------------------------------
*/
static int mouse_down[3];
static int omx, omy, mx, my;
static void get_from_UI()
{

    if (!mouse_down[0] && !mouse_down[2]) return;

    float i = ((mx / (float)win_x) * NUM_COL);
    float j = (((win_y - my) / (float)win_y) * NUM_ROW);

    if (i<1 || i>NUM_COL || j<1 || j>NUM_ROW) return;

    if (mouse_down[0]) {
        f3 f = { fluidAmount * (mx - omx), fluidAmount * (omy - my),  fluidAmount * (mx - omx + omy - my)/2 };
        f3 pos = { i, j };//{ NUM_COL / 2, NUM_ROW / 2 }; //
        addInk(f, pos, fluidSize);
    }

    if (mouse_down[2]) {
        f3 f = { forceAmount * (mx - omx), forceAmount * (omy - my) };
        f3 pos = { i, j };//{ NUM_COL / 2, NUM_ROW / 2 }; //
        addExternalForce(f, pos, forceSize);
    }

    omx = mx;
    omy = my;

    return;
}

static void mouse_func(int button, int state, int x, int y)
{
    omx = mx = x;
    omx = my = y;

    mouse_down[button] = state == GLUT_DOWN;
}

static void motion_func(int x, int y)
{
    mx = x;
    my = y;
}

static void key_func(unsigned char key, int x, int y)
{
    switch (key)
    {
    //case 'c':
    //case 'C':
    //    clear_data();
    //    break;

    //case 'q':
    //case 'Q':
    //    free_data();
    //    exit(0);
    //    break;

    case 'v':
    case 'V':
        dvel = !dvel;
        break;

    //case 'r':
    //case 'R':
    //    N = 128;
    //    dt = 0.1f;
    //    diff = 0.0f;
    //    visc = 0.0f;
    //    force = 0.5f;
    //    source = 100.0f;
    //    break;

    case 'a':
        dt += dt_del;
        printf("dt is now %f\n", dt);
        break;

    case 'A':
        dt -= dt_del;
        printf("dt is now %f\n", dt);
        break;

    case 's':
        visc += visc_del;
        printf("visc is now %f\n", visc);
        break;

    case 'S':
        visc -= visc_del;
        printf("visc is now %f\n", visc);
        break;

    case 'd':
        vort += vort_del;
        printf("vort is now %f\n", vort);
        break;

    case 'D':
        vort -= vort_del;
        printf("vort is now %f\n", vort);
        break;

    case 'f':
        fluidSize += fluidSize_del;
        printf("fluidSize is now %f\n", fluidSize);
        break;

    case 'F':
        fluidSize -= fluidSize_del;
        printf("fluidSize is now %f\n", fluidSize);
        break;

    case 'g':
        fluidAmount += fluidAmount_del;
        printf("fluidAmount is now %f\n", fluidAmount);
        break;

    case 'G':
        fluidAmount -= fluidAmount_del;
        printf("fluidAmount is now %f\n", fluidAmount);
        break;

    case 'h':
        forceSize += forceSize_del;
        printf("forceSize is now %f\n", forceSize);
        break;

    case 'H':
        forceSize -= forceSize_del;
        printf("forceSize is now %f\n", forceSize);
        break;

    case 'j':
        forceAmount += forceAmount_del;
        printf("forceAmount is now %f\n", forceAmount);
        break;

    case 'J':
        forceAmount -= forceAmount_del;
        printf("forceAmount is now %f\n", forceAmount);
        break;

        // case 'p':
        // case 'P':
        // 	if (square == 0)
        // 	{
        // 		square = 1;
        // 	}
        // 	else
        // 	{
        // 		square = 0;
        // 	}
        // 	break;
    }

    // Stay within bounds.
    if (dt > dtMax)
        dt = dtMax;
    if (dt < dtMin)
        dt = dtMin;
    if (visc > viscMax)
        visc = viscMax;
    if (visc < viscMin)
        visc = viscMin;
    if (vort > vortMax)
        vort = vortMax;
    if (vort < vortMin)
        vort = vortMin;
    if (fluidSize > fluidSizeMax)
        fluidSize = fluidSizeMax;
    if (fluidSize < fluidSizeMin)
        fluidSize = fluidSizeMin;
    if (fluidAmount > fluidAmountMax)
        fluidAmount = fluidAmountMax;
    if (fluidAmount < fluidAmountMin)
        fluidAmount = fluidAmountMin;
    if (forceSize > forceSizeMax)
        forceSize = forceSizeMax;
    if (forceSize < forceSizeMin)
        forceSize = forceSizeMin;
    if (forceAmount > forceAmountMax)
        forceAmount = forceAmountMax;
    if (forceAmount < forceAmountMin)
        forceAmount = forceAmountMin;
}

/*
  ----------------------------------------------------------------------
   open_glut_window --- open a glut compatible window and set callbacks
  ----------------------------------------------------------------------
*/

static void open_glut_window(void)
{
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

	glutInitWindowPosition(0, 0);
	glutInitWindowSize(win_x, win_y);
	win_id = glutCreateWindow("Fluid Simulation");

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();

	pre_display();
    glutKeyboardFunc(key_func);
    glutMouseFunc(mouse_func);
    glutMotionFunc(motion_func);
	glutIdleFunc(idle_func);
	glutDisplayFunc(display_func);
}


/*
  ----------------------------------------------------------------------
   main --- main routine
  ----------------------------------------------------------------------
*/

int main(int argc, char** argv)
{
	glutInit(&argc, argv);

	if (argc != 1 && argc != 6) {
		fprintf(stderr, "usage : %s N dt diff visc force source\n", argv[0]);
		fprintf(stderr, "where:\n"); \
			fprintf(stderr, "\t N      : grid resolution\n");
		fprintf(stderr, "\t dt     : time step\n");
		fprintf(stderr, "\t diff   : diffusion rate of the density\n");
		fprintf(stderr, "\t visc   : viscosity of the fluid\n");
		fprintf(stderr, "\t force  : scales the mouse movement that generate a force\n");
		fprintf(stderr, "\t source : amount of density that will be deposited\n");
		exit(1);
	}

	if (argc == 1) {
		dt = 0.1f;
		diff = 0.0f;
		visc = 0.0f;
		fprintf(stderr, "Using defaults : dt=%g diff=%g visc=%g\n",
			 dt, diff, visc);
	}
	else {
		dt = atof(argv[2]);
		diff = atof(argv[3]);
		visc = atof(argv[4]);
	}

    // Setting added vars
    dt = 0.1;
    visc = 0.01;
    vort = 0.5;
    fluidSize = 0.5;
    fluidAmount = 0.5;
    forceSize = 0.5;
    forceAmount = 0.5;
    //

    // simulate initial hit, update u by inrementing the timestep
    

    clearGrids();

	win_x = 512;
	win_y = 512;
	open_glut_window();

	glutMainLoop();

	exit(0);
}


void update_grid(f3 u[NUM_ROW][NUM_COL], f3 u_temp[NUM_ROW][NUM_COL], f3 p[NUM_ROW][NUM_COL], f3 divergence_u[NUM_ROW][NUM_COL], float timestep, float epsilon) {

<<<<<<< HEAD
    // Apply the first 3 operators in Equation 12. 
    advect(u_temp, u, u, timestep, 1);

    advect(ink[(q_ink + 1) % 1], u, ink[q_ink % 1], timestep, 1);
    q_ink = (q_ink + 1) % 1;

    //u1 will store the output of diffuse
    diffuse(u, u_temp, epsilon, timestep);

    addForces(u, timestep);

    //for (int i = 0; i < NUM_ROW; i++) {
    //    for (int j = 0; j < NUM_COL; j++) {
    //        u_temp[j][i] = u[j][i];
    //    }
    //}

    divergence(divergence_u, u);

 
    computePressure(p, u, divergence_u, timestep);

    subtractPressureGradient(u_temp, u, p);
=======
    // velocity update
    addForces(u, timestep);
    diffuse(u_temp, u, epsilon, timestep);
    divergence(divergence_u, u);
    computePressure(p, u, divergence_u, timestep);
    subtractPressureGradient(u_temp, u, p);
    advect(u, u_temp, u_temp, timestep, 1);

    // density step
    diffuse(ink[(q_ink + 1) % 1], ink[q_ink % 1], epsilon, timestep);
    q_ink = (q_ink + 1) % 1;

    advect_density(ink[(q_ink + 1) % 1], u_temp, ink[q_ink % 1], timestep, 1);
    q_ink = (q_ink + 1) % 1;
>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62
}

//void advect(float2 coords : WPOS,   // grid coordinates     
//    out float4 xNew : COLOR,  // advected qty     
//    uniform float timestep,             
//    uniform    float rdx,        // 1 / grid scale    
//    uniform    samplerRECT u,    // input velocity     
//    uniform    samplerRECT x)    // qty to advect 
//{   // follow the velocity field "back in time"     
//    float2 pos = coords - timestep * rdx * f2texRECT(u, coords); 
// // interpolate and write to the output fragment   
//    xNew = f4texRECTbilerp(x, pos); 
//} 
void advect(f3 u[NUM_ROW][NUM_COL], f3 u1[NUM_ROW][NUM_COL], f3 X[NUM_ROW][NUM_COL], float timestep, float gain) {

    // inner elements
    for (int x = 1; x < NUM_COL - 1; x++) {
        for (int y = 1; y < NUM_ROW - 1; y++) {
            // follow the velocity field "back in time" 
            float x_old = (x - (u1[y][x].x * timestep * NUM_COL)); if (x_old < 0.5) x_old = 0.5; if (x_old > NUM_COL - 1.5) x_old = NUM_COL - 1.5;
            float y_old = (y - (u1[y][x].y * timestep * NUM_ROW)); if (y_old < 0.5) y_old = 0.5; if (y_old > NUM_ROW - 1.5) y_old = NUM_ROW - 1.5;

            
            int y_1 = (int)(y_old);
            int y_2 = y_1 + 1;
            int x_1 = (int)(x_old);
            int x_2 = x_1 + 1;


            // interpolate and update u
            float p_x = x_old - (int)x_old, p_y = y_old - (int)y_old;

            u[y][x] = (1 - p_y) * ((1 - p_x) * X[y_1][x_1] + (p_x)*X[y_1][x_2]) + (p_y) * ((1 - p_x) * X[y_2][x_1] + (p_x)*X[y_2][x_2]);
            u[y][x] = u[y][x] * gain;

        }
    }

    // upper + lower border
    for (int x = 0; x < NUM_COL; x++) {
        //upper border
<<<<<<< HEAD
        u[0][x] =  -u[1][x];

        //lower border
        u[NUM_ROW - 1][x] = -u[NUM_ROW - 2][x];
=======
        u[0][x].x = -u[1][x].x;
        u[0][x].y = u[1][x].y;

        //lower border
        u[NUM_ROW - 1][x].x = -u[NUM_ROW - 2][x].x;
        u[NUM_ROW - 1][x].y = u[NUM_ROW - 2][x].y;
>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62
    }

    //left + right border
    for (int y = 0; y < NUM_ROW; y++) {
        //left border
<<<<<<< HEAD
        u[y][0] = -u[y][1];

        //right border
        u[y][NUM_COL - 1] = -u[y][NUM_COL - 2];
    }

    // corners
    u[0][0]                     = (u[0][1] + u[1][0]) / 3;
    u[NUM_ROW - 1][0]           = (u[NUM_ROW - 1][1] + u[NUM_ROW - 2][0]) / 3;
    u[0][NUM_COL - 1]           = (u[0][NUM_COL - 2] + u[1][NUM_COL - 1]) / 3;
    u[NUM_ROW - 1][NUM_COL - 1]  = (u[NUM_ROW - 1][NUM_COL - 2] + u[NUM_ROW - 2][NUM_COL - 1]) / 3;
=======
        u[y][0].x = u[y][1].x;
        u[y][0].y = -u[y][1].y;

        //right border
        u[y][NUM_COL - 1].x = u[y][NUM_COL - 2].x;
        u[y][NUM_COL - 1].y = -u[y][NUM_COL - 2].y;
    }

    u[0][0] = 0.5f * (u[0][1] + u[1][0]);
    u[NUM_ROW - 1][0] = 0.5f * (u[NUM_ROW - 1][1] + u[NUM_ROW - 2][0]);
    u[0][NUM_COL - 1] = 0.5f * (u[0][NUM_COL - 2] + u[1][NUM_COL - 1]);
    u[NUM_ROW - 1][NUM_COL - 1] = 0.5f * (u[NUM_ROW - 1][NUM_COL - 2] + u[NUM_ROW - 2][NUM_COL - 1]);
}

void advect_density(f3 u[NUM_ROW][NUM_COL], f3 u1[NUM_ROW][NUM_COL], f3 X[NUM_ROW][NUM_COL], float timestep, float gain) {

    // inner elements
    for (int x = 1; x < NUM_COL - 1; x++) {
        for (int y = 1; y < NUM_ROW - 1; y++) {
            // follow the velocity field "back in time" 
            float x_old = (x - (u1[y][x].x * timestep * NUM_COL)); if (x_old < 0.5) x_old = 0.5; if (x_old > NUM_COL - 1.5) x_old = NUM_COL - 1.5;
            float y_old = (y - (u1[y][x].y * timestep * NUM_ROW)); if (y_old < 0.5) y_old = 0.5; if (y_old > NUM_ROW - 1.5) y_old = NUM_ROW - 1.5;


            int y_1 = (int)(y_old);
            int y_2 = y_1 + 1;
            int x_1 = (int)(x_old);
            int x_2 = x_1 + 1;


            // interpolate and update u
            float p_x = x_old - (int)x_old, p_y = y_old - (int)y_old;

            u[y][x] = (1 - p_y) * ((1 - p_x) * X[y_1][x_1] + (p_x)*X[y_1][x_2]) + (p_y) * ((1 - p_x) * X[y_2][x_1] + (p_x)*X[y_2][x_2]);
        }
    }

    // upper + lower border
    for (int x = 0; x < NUM_COL; x++) {
        //upper border
        u[0][x].x = u[1][x].x;
        u[0][x].y = u[1][x].y;

        //lower border
        u[NUM_ROW - 1][x].x = u[NUM_ROW - 2][x].x;
        u[NUM_ROW - 1][x].y = u[NUM_ROW - 2][x].y;
    }

    //left + right border
    for (int y = 0; y < NUM_ROW; y++) {
        //left border
        u[y][0].x = u[y][1].x;
        u[y][0].y = u[y][1].y;

        //right border
        u[y][NUM_COL - 1].x = u[y][NUM_COL - 2].x;
        u[y][NUM_COL - 1].y = u[y][NUM_COL - 2].y;
    }

    u[0][0] = 0.5f * (u[0][1] + u[1][0]);
    u[NUM_ROW - 1][0] = 0.5f * (u[NUM_ROW - 1][1] + u[NUM_ROW - 2][0]);
    u[0][NUM_COL - 1] = 0.5f * (u[0][NUM_COL - 2] + u[1][NUM_COL - 1]);
    u[NUM_ROW - 1][NUM_COL - 1] = 0.5f * (u[NUM_ROW - 1][NUM_COL - 2] + u[NUM_ROW - 2][NUM_COL - 1]);
>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62
}

void jacobi_borderConditions(f3 u[NUM_ROW][NUM_COL], int sign) {
    // upper + lower border
    for (int x = 0; x < NUM_COL; x++) {
        //upper border
        u[0][x] = -u[1][x] * sign;

        //lower border
        u[NUM_ROW - 1][x] = -u[NUM_ROW - 2][x] * sign;
    }

    //left + right border
    for (int y = 0; y < NUM_ROW; y++) {
        //left border
        u[y][0] = -u[y][1] * sign;

        //right border
        u[y][NUM_COL - 1] = -u[y][NUM_COL - 2] * sign;
    }
}

void jacobi_diffuse(f3 xNew[NUM_ROW][NUM_COL], f3 guess[NUM_ROW][NUM_COL], f3 b[NUM_ROW][NUM_COL], int numIterations, float dt) {

    for (int x = 1; x < NUM_COL - 1; x++) {
        for (int y = 1; y < NUM_ROW -1; y++) {
            f3 alpha = (guess[y][x] * guess[y][x] / dt);
            xNew[y][x] = (guess[y - 1][x] + guess[y + 1][x] + guess[y][x - 1] + guess[y][x + 1] + alpha * b[y][x]) * 1 / (4 + alpha);
        }
        jacobi_borderConditions(xNew, 1);
    }
 

    int iter;
    for (iter = 0; iter < numIterations; iter++) {
        if (iter % 2 == 0) {
            for (int x = 1; x < NUM_COL - 1; x++) {
                for (int y = 1; y < NUM_ROW - 1; y++) {
                    f3 alpha = (xNew[y][x] * xNew[y][x] / dt);
                    guess[y][x] = (xNew[y - 1][x] + xNew[y + 1][x] + xNew[y][x - 1] + xNew[y][x + 1] + alpha * b[y][x]) * 1 / (4 + alpha);
                }
            }
            jacobi_borderConditions(guess, 1);
        }
        else {
            for (int x = 1; x < NUM_COL - 1; x++) {
                for (int y = 1; y < NUM_ROW - 1; y++) {
                    f3 alpha = (guess[y][x] * guess[y][x] / dt);
                    xNew[y][x] = (guess[y - 1][x] + guess[y + 1][x] + guess[y][x - 1] + guess[y][x + 1] + alpha * b[y][x]) * 1 / (4 + alpha);
                }
            }
            jacobi_borderConditions(xNew, 1);
        }

    }

    //need to copy result to xNew
    if (iter % 2 == 0) {
        for (int x = 1; x < NUM_COL - 1; x++) {
            for (int y = 1; y < NUM_ROW - 1; y++) {
                xNew[y][x] = guess[y][x];
            }
        }
    }
}

void diffuse(f3 u[NUM_ROW][NUM_COL], f3 u1[NUM_ROW][NUM_COL], float diffusion_rate, float timestep) {

<<<<<<< HEAD
    //compute inner elements
    jacobi_diffuse(u, u1, u1, 50, timestep * diffusion_rate);
=======
    for (int x = 0; x < NUM_COL; x++) {
        for (int y = 0; y < NUM_ROW; y++) {
            U_b[y][x] = u1[y][x];
        }
    }

    //compute inner elements
    jacobi_diffuse(u, u1, U_b, diffusion_iter, timestep * diffusion_rate);
>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62
}

void addForces(f3 u[NUM_ROW][NUM_COL], float timestep) {
    GUI_force* f;
    f3 impulse;
    f3 tmp;
    while (!queueEmpty()) {
        f = force_dequeue();
        impulse = f->force * timestep;
<<<<<<< HEAD
=======
        int x = (float)f->position.x;
        int y = (float)(f->position.y);
        //u[y][x] = u[y][x] + impulse;
>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62
        for (int x = 1; x < NUM_COL - 1; x++) {
            for (int y = 1; y < NUM_ROW - 1; y++) {
                f3 pt = { (float)x, (float)y };
                
                float effect = (float)exp(-pt.dist(f->position) / f->impulseRadius);
                tmp = (effect * impulse);
                u[y][x] = u[y][x] + tmp;
            }
        }
        free(f);
    }

    int inkUpdated = 0;
    while (!ink_queueEmpty()) {
        inkUpdated = 1;
        f = ink_dequeue();
        impulse = f->force * timestep;
        for (int x = 1; x < NUM_COL - 1; x++) {
            for (int y = 1; y < NUM_ROW - 1; y++) {
                f3 pt = { (float)x, (float)y };
                float effect = (float)exp(-pt.dist(f->position) / f->impulseRadius);
                tmp = (effect * impulse);
                ink[(q_ink + 1) % 1][y][x] = ink[(q_ink + 1) % 1][y][x] + f3::abs(tmp);
<<<<<<< HEAD
=======
                ink[(q_ink + 1) % 1][y][x].cap(1);
>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62
            }
        }
        free(f);
    }
    if (inkUpdated) {
        q_ink = (q_ink + 1) % 1;
    }
}

void divergence(f3 div_f[NUM_ROW][NUM_COL], f3 f[NUM_ROW][NUM_COL]) {
    // inner elements
    for (int x = 1; x < NUM_COL; x++) {
        for (int y = 1; y < NUM_ROW; y++) {
<<<<<<< HEAD
            div_f[y][x].x = (f[y][x + 1].x - f[y][x - 1].x) / (float)2.0;
            div_f[y][x].x = (f[y - 1][x].y - f[y + 1][x].y) / (float)2.0;
=======
            div_f[y][x].x = (f[y][x + 1].x - f[y][x - 1].x) / (float)(-2.0 * NUM_COL);
            div_f[y][x].y = (f[y + 1][x].y - f[y - 1][x].y) / (float)(-2.0 * NUM_ROW);
>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62
        }
    }

    // upper + lower border
    for (int x = 0; x < NUM_COL; x++) {
        //upper border
        div_f[0][x].x = (f[1][x].x) / (float)2.0;
        div_f[0][x].y = (-f[1][x].y) / (float)2.0;

        //lower border
        div_f[NUM_ROW - 1][x].x = (-f[NUM_ROW - 2][x].x) / (float)2.0;
        div_f[NUM_ROW - 1][x].y = (f[NUM_ROW - 2][x].y) / (float)2.0;
    }

    //left + right border
    for (int y = 0; y < NUM_ROW; y++) {
        //left border
        div_f[y][0].x = (f[y][1].x) / (float)2.0;
        div_f[y][0].y = (-f[y][1].y) / (float)2.0;

        //right border
        div_f[y][NUM_COL - 1].x = (-f[y][NUM_COL - 2].x) / (float)2.0;
        div_f[y][NUM_COL - 1].y = (f[y][NUM_COL - 2].y) / (float)2.0;
    }
}

void jacobi_pressure(f3 xNew[NUM_ROW][NUM_COL], f3 X[NUM_ROW][NUM_COL], f3 b[NUM_ROW][NUM_COL], int numIterations, float dt) {
    for (int x = 0; x < NUM_COL; x++) {
        for (int y = 0; y < NUM_ROW; y++) {
            xNew[y][x] = f3(0, 0);
        }
    }

    int iter;
    for (iter = 1; iter < numIterations; iter++) {
        if (iter % 2 == 0) {
            for (int x = 1; x < NUM_COL - 1; x++) {
                for (int y = 1; y < NUM_ROW - 1; y++) {
                    f3 alpha = -xNew[y][x] * xNew[y][x];
                    X[y][x] = (xNew[y - 1][x] + xNew[y + 1][x] + xNew[y][x - 1] + xNew[y][x + 1] + alpha * b[y][x]) * 0.25;
                }
            }
            jacobi_borderConditions(X, -1);
        }
        else {
            for (int x = 1; x < NUM_COL - 1; x++) {
                for (int y = 1; y < NUM_ROW - 1; y++) {
                    f3 alpha = -X[y][x] * X[y][x];
                    xNew[y][x] = (X[y - 1][x] + X[y + 1][x] + X[y][x - 1] + X[y][x + 1] + alpha * b[y][x]) * 0.25;
                }
            }
            jacobi_borderConditions(xNew, -1);
        }
    }

    //need to copy result to xNew
    if (iter % 2 == 0) {
        for (int x = 1; x < NUM_COL - 1; x++) {
            for (int y = 1; y < NUM_ROW - 1; y++) {
                xNew[y][x] = X[y][x];
            }
        }
    }
}


void computePressure(f3 p[NUM_ROW][NUM_COL], f3 u[NUM_ROW][NUM_COL], f3 divergence_u[NUM_ROW][NUM_COL], float timestep) {
<<<<<<< HEAD
    jacobi_pressure(p, P_guess, divergence_u, 40, timestep);

    //// upper + lower border
    //for (int x = 0; x < NUM_COL; x++) {
    //    //upper border
    //    p[0][x] = p[1][x];

    //    //lower border
    //    p[NUM_ROW - 1][x] = p[NUM_ROW - 2][x];
    //}

    ////left + right border
    //for (int y = 0; y < NUM_ROW; y++) {
    //    //left border
    //    p[y][0] = p[y][1];

    //    //right border
    //    p[y][NUM_COL - 1] = p[y][NUM_COL - 2];
    //}
=======
    jacobi_pressure(p, P_guess, divergence_u, pressure_iter, timestep);


>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62
}


void subtractPressureGradient(f3 u[NUM_ROW][NUM_COL], f3 u1[NUM_ROW][NUM_COL], f3 p[NUM_ROW][NUM_COL]) {
    for (int y = 1; y < NUM_ROW - 1; y++) {
        for (int x = 1; x < NUM_COL - 1; x++) {
<<<<<<< HEAD
            u[y][x] = u1[y][x] - f3(p[y][x + 1].x - p[y][x - 1].x, p[y - 1][x].y - p[y + 1][x].y);
        }
    }
=======
            u[y][x] = u1[y][x] - 0.5 * NUM_COL * f3(p[y][x + 1].x - p[y][x - 1].x, p[y + 1][x].y - p[y - 1][x].y);
        }
    }

    // upper + lower border
    for (int x = 0; x < NUM_COL; x++) {
        //upper border
        u[0][x].x = -u[1][x].x;
        u[0][x].y = u[1][x].y;

        //lower border
        u[NUM_ROW - 1][x].x = -u[NUM_ROW - 2][x].x;
        u[NUM_ROW - 1][x].y = u[NUM_ROW - 2][x].y;
    }

    //left + right border
    for (int y = 0; y < NUM_ROW; y++) {
        //left border
        u[y][0].x = u[y][1].x;
        u[y][0].y = -u[y][1].y;

        //right border
        u[y][NUM_COL - 1].x = u[y][NUM_COL - 2].x;
        u[y][NUM_COL - 1].y = -u[y][NUM_COL - 2].y;
    }

    u[0][0] = 0.5f * (u[0][1] + u[1][0]);
    u[NUM_ROW - 1][0] = 0.5f * (u[NUM_ROW - 1][1] + u[NUM_ROW - 2][0]);
    u[0][NUM_COL - 1] = 0.5f * (u[0][NUM_COL - 2] + u[1][NUM_COL - 1]);
    u[NUM_ROW - 1][NUM_COL - 1] = 0.5f * (u[NUM_ROW - 1][NUM_COL - 2] + u[NUM_ROW - 2][NUM_COL - 1]);
>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62
}

/**
* TODO: Should add forces to a queue
*/
void addExternalForce(f3 force, f3 position, float impulseRadius) {
    force_enqueue(force, position, impulseRadius);
}

void addInk(f3 force, f3 position, float impulseRadius) {
    ink_enqueue(force, position, impulseRadius);
<<<<<<< HEAD
}
=======
}
>>>>>>> 2b31b9ee5a02e628a1cd5ebae7c2189ae02d3a62
