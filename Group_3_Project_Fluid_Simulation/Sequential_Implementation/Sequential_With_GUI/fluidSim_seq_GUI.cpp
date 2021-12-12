#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <chrono>
#include <iostream>

#include "freeglut.h"
#include "utilities.h"


#define NUM_ROW 128
#define NUM_COL NUM_ROW
#define GRID_SIZE NUM_ROW * NUM_COL
#define BOUNDARY_GAIN 0.75

static void get_from_UI();
static void display_func(void);

// update functions
void forces(f3 u[NUM_ROW][NUM_COL], float timestep);
void advect(f3 f[NUM_ROW][NUM_COL], f3 x[NUM_ROW][NUM_COL], f3 f_prev[NUM_ROW][NUM_COL], float timestep, float gain, f3 borderConditions);
void diffuse(f3 x[NUM_ROW][NUM_COL], f3 x_prev[NUM_ROW][NUM_COL], float diffusion_rate, float dt, f3 borderConditions);
void pressure(float p[NUM_ROW][NUM_COL], f3 u[NUM_ROW][NUM_COL], float timestep);
void setBorders(f3 x[NUM_ROW][NUM_COL], int borderCondition);
void setBorders(float x[NUM_ROW][NUM_COL], int borderCondition);
void update_grid(float timestep, float epsilon);

// timing variables
auto begin = std::chrono::system_clock::now();
auto end = std::chrono::system_clock::now();


// GUI event
void addExternalForce(f3 force, f3 position, float impulseRadius);
void addInk(f3 force, f3 position, float impulseRadius);


// assign pointers for u, u1, u2 2D grids
f3 U[NUM_ROW][NUM_COL], U_prev[NUM_ROW][NUM_COL]; //velocity fields
float P[NUM_ROW][NUM_COL], P_prev[NUM_ROW][NUM_COL]; //pressure fields
f3 ink[NUM_ROW][NUM_COL], ink_prev[NUM_ROW][NUM_COL]; //density fields

volatile int dvel, show_commands = 1, shoot_liquid = 0; // used to toggle displays between velocity and density

static int win_id;
static int win_x, win_y;
volatile static float ;

static int diffusion_iter = 100, pressure_iter = 50;

volatile static float dt, diff, visc, vort, fluidSize, fluidAmount, forceSize, forceAmount;

// Variable bounds.
float numIntervals = 100;
float dtMax = 5;
float dtMin = 0.01;
float dt_del = (dtMax - dtMin) / numIntervals;
float viscMax = 1;
float viscMin = 1e-4;
float visc_del = (viscMax - viscMin) / numIntervals;
float vortMax = 1.;
float vortMin = 1e-6;
float vort_del = (vortMax - vortMin) / numIntervals;
float fluidSizeMax = 5.;
float fluidSizeMin = 0.1;
float fluidSize_del = (fluidSizeMax - fluidSizeMin) / numIntervals;
float fluidAmountMax = 3.;
float fluidAmountMin = 0.01;
float fluidAmount_del = (fluidAmountMax - fluidAmountMin) / numIntervals;
float forceSizeMax = 5.;
float forceSizeMin = 0.1;
float forceSize_del = (forceSizeMax - forceSizeMin) / numIntervals;
float forceAmountMax = 40.;
float forceAmountMin = 5;
float forceAmount_del = (forceAmountMax - forceAmountMin) / numIntervals;





void clearGrids() {
    for (int x = 0; x < NUM_COL; x++) {
        for (int y = 0; y < NUM_ROW; y++) {
            P[y][x] = P_prev[y][x] = 0;
            U[y][x] = U_prev[y][x] = ink[y][x] = ink_prev[y][x] = f3(0,0);
        }
    }
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


// draw velocity field
static void draw_velocity(void)
{
	int i, j;
	float x, y, h;

	h = 1.0f / NUM_COL;

	glColor3f(1.0f, 1.0f, 1.0f);
	glLineWidth(1.0f);

	glBegin(GL_LINES);

	for (i = 0; i < NUM_COL; i++) {
		x = (i - 0.5f) * h;
		for (j = 0; j < NUM_ROW; j++) {
			y = (j - 0.5f) * h;

			glVertex2f(x, y);
			glVertex2f(x + U[j][i].x, y + U[j][i].y);
		}
	}

	glEnd();
}

// draw density field
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

            d00 = ink[j][i];
            d01 = ink[j + 1][i];
            d10 = ink[j][i + 1];
            d11 = ink[j + 1][i + 1];

            glColor3f(d00.x, d00.y, d00.z); glVertex2f(x, y);
            glColor3f(d10.x, d10.y, d10.z); glVertex2f(x + h, y);
            glColor3f(d11.x, d11.y, d11.z); glVertex2f(x + h, y + h);
            glColor3f(d01.x, d01.y, d01.z); glVertex2f(x, y + h);

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
    update_grid(dt, visc);

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

    //get framerate
    end = std::chrono::system_clock::now();
    auto dur = end - begin;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    printf("framerate: %f\n", 1000.0 / ms);
    begin = std::chrono::system_clock::now();

    //display sliders
    if (show_commands) Sliders();

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
    if (shoot_liquid) {
        float i = (NUM_ROW) / 5;
        float j = (NUM_COL) / 5;

        f3 f = { forceAmount, forceAmount};
        f3 pos = { i, j };
        addExternalForce(f, pos, forceSize);
        
        f = { fluidAmount, fluidAmount,  0.0f};
        pos = { i, j };
        addInk(f, pos, fluidSize);
    }

    if (!mouse_down[0] && !mouse_down[2]) return;

    float i = ((mx / (float)win_x) * NUM_COL);
    float j = (((win_y - my) / (float)win_y) * NUM_ROW);

    if (i<1 || i>NUM_COL || j<1 || j>NUM_ROW) return;

    if (mouse_down[0]) {
        f3 f = { fluidAmount * (mx - omx), fluidAmount * (omy - my),  fluidAmount * (mx - omx + omy - my)/2 };
        f3 pos = { i, j };
        addInk(f, pos, fluidSize);
    }

    if (mouse_down[2]) {
        f3 f = { forceAmount * (mx - omx), forceAmount * (omy - my) };
        f3 pos = { i, j };
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
    case 'c':
    case 'C':
        clearGrids();
        break;

    case 'q':
    case 'Q':
        show_commands = !show_commands;
        break;
    case 'w':
    case 'W':
        shoot_liquid = !shoot_liquid;
        break;
    case 'v':
    case 'V':
        dvel = !dvel;
        break;

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
    fluidSize = 3;
    fluidAmount = 3;
    forceSize = 2;
    forceAmount = 40;
    //

    // simulate initial hit, update u by inrementing the timestep
    

    clearGrids();

	win_x = 800;
	win_y = 800;
	open_glut_window();

	glutMainLoop();

	exit(0);
}


void update_grid(float timestep, float epsilon) {

    forces(U, timestep);

    // velocity update
    diffuse(U_prev, U, 0.0, timestep, f3(1, 2));
    pressure(P, U_prev, timestep);
    advect(U, U_prev, U_prev, timestep, 1, f3(1, 2));
    pressure(P, U, timestep);
    
    // density step
    diffuse(ink_prev, ink, 0.0, timestep, f3(3, 3)); //TODO: have dynamic memory allocation
    advect(ink, U, ink_prev, timestep, 1, f3(3, 3));
}



void setBorders(f3 x[NUM_ROW][NUM_COL], int borderCondition) {
    // upper + lower border
    for (int ixd = 0; ixd < NUM_COL; ixd++) {
        switch (borderCondition) {
        case 1:
            x[0][ixd]               = -x[1][ixd]; //upper border
            x[NUM_ROW - 1][ixd]     = -x[NUM_ROW - 2][ixd]; //lower border
            x[ixd][0]               = x[ixd][1]; //left border
            x[ixd][NUM_COL - 1]     = x[ixd][NUM_COL - 2]; //right border
            break;
        case 2:
            x[0][ixd]               = x[1][ixd]; //upper border
            x[NUM_ROW - 1][ixd]     = x[NUM_ROW - 2][ixd]; //lower border
            x[ixd][0]               = -x[ixd][1]; //left border
            x[ixd][NUM_COL - 1]     = -x[ixd][NUM_COL - 2]; //right border
            break;
        default:
            x[0][ixd]               = x[1][ixd]; //upper border
            x[NUM_ROW - 1][ixd]     = x[NUM_ROW - 2][ixd]; //lower border
            x[ixd][0]               = x[ixd][1]; //left border
            x[ixd][NUM_COL - 1]     = x[ixd][NUM_COL - 2]; //right border
            break;
        }
        
    }

    x[0][0] = 0.5f * (x[0][1] + x[1][0]);
    x[NUM_ROW - 1][0] = 0.5f * (x[NUM_ROW - 1][1] + x[NUM_ROW - 2][0]);
    x[0][NUM_COL - 1] = 0.5f * (x[0][NUM_COL - 2] + x[1][NUM_COL - 1]);
    x[NUM_ROW - 1][NUM_COL - 1] = 0.5f * (x[NUM_ROW - 1][NUM_COL - 2] + x[NUM_ROW - 2][NUM_COL - 1]);
}

void setBorders(f3 x[NUM_ROW][NUM_COL], int borderConditionX, int borderConditionY) {
    // upper + lower border
    for (int ixd = 0; ixd < NUM_COL; ixd++) {
        switch (borderConditionX) {
        case 1:
            x[0][ixd].x = -x[1][ixd].x;  //upper border
            x[NUM_ROW - 1][ixd].x = -x[NUM_ROW - 2][ixd].x;  //lower border
            x[ixd][0].x = x[ixd][1].x;  //left border
            x[ixd][NUM_COL - 1].x = x[ixd][NUM_COL - 2].x;  //right border
            break;
        case 2:
            x[0][ixd].x = x[1][ixd].x;  //upper border
            x[NUM_ROW - 1][ixd].x = x[NUM_ROW - 2][ixd].x;  //lower border
            x[ixd][0].x = -x[ixd][1].x;  //left border
            x[ixd][NUM_COL - 1].x = -x[ixd][NUM_COL - 2].x;  //right border
            break;
        default:
            x[0][ixd].x = x[1][ixd].x;  //upper border
            x[NUM_ROW - 1][ixd].x = x[NUM_ROW - 2][ixd].x;  //lower border
            x[ixd][0].x = x[ixd][1].x;  //left border
            x[ixd][NUM_COL - 1].x = x[ixd][NUM_COL - 2].x;  //right border
            break;
        }
        switch (borderConditionY) {
        case 1:
            x[0][ixd].y  = -x[1][ixd].y;  //upper border
            x[NUM_ROW - 1][ixd].y  = -x[NUM_ROW - 2][ixd].y;  //lower border
            x[ixd][0].y  = x[ixd][1].y;  //left border
            x[ixd][NUM_COL - 1].y  = x[ixd][NUM_COL - 2].y;  //right border
            break;
        case 2:
            x[0][ixd].y  = x[1][ixd].y;  //upper border
            x[NUM_ROW - 1][ixd].y  = x[NUM_ROW - 2][ixd].y;  //lower border
            x[ixd][0].y  = -x[ixd][1].y;  //left border
            x[ixd][NUM_COL - 1].y  = -x[ixd][NUM_COL - 2].y;  //right border
            break;
        default:
            x[0][ixd].y  = x[1][ixd].y;  //upper border
            x[NUM_ROW - 1][ixd].y  = x[NUM_ROW - 2][ixd].y;  //lower border
            x[ixd][0].y  = x[ixd][1].y;  //left border
            x[ixd][NUM_COL - 1].y  = x[ixd][NUM_COL - 2].y;  //right border
            break;
        }
    }

    x[0][0] = 0.5f * (x[0][1] + x[1][0]);
    x[NUM_ROW - 1][0] = 0.5f * (x[NUM_ROW - 1][1] + x[NUM_ROW - 2][0]);
    x[0][NUM_COL - 1] = 0.5f * (x[0][NUM_COL - 2] + x[1][NUM_COL - 1]);
    x[NUM_ROW - 1][NUM_COL - 1] = 0.5f * (x[NUM_ROW - 1][NUM_COL - 2] + x[NUM_ROW - 2][NUM_COL - 1]);
}

void setBorders(float x[NUM_ROW][NUM_COL], int borderCondition) {
    // upper + lower border
    for (int ixd = 0; ixd < NUM_COL; ixd++) {
        switch (borderCondition) {
        case 1:
            x[0][ixd] = -x[1][ixd]; //upper border
            x[NUM_ROW - 1][ixd] = -x[NUM_ROW - 2][ixd]; //lower border
            x[ixd][0] = x[ixd][1]; //left border
            x[ixd][NUM_COL - 1] = x[ixd][NUM_COL - 2]; //right border
            break;
        case 2:
            x[0][ixd] = x[1][ixd]; //upper border
            x[NUM_ROW - 1][ixd] = x[NUM_ROW - 2][ixd]; //lower border
            x[ixd][0] = -x[ixd][1]; //left border
            x[ixd][NUM_COL - 1] = -x[ixd][NUM_COL - 2]; //right border
            break;
        default:
            x[0][ixd] = x[1][ixd]; //upper border
            x[NUM_ROW - 1][ixd] = x[NUM_ROW - 2][ixd]; //lower border
            x[ixd][0] = x[ixd][1]; //left border
            x[ixd][NUM_COL - 1] = x[ixd][NUM_COL - 2]; //right border
            break;
        }

    }

    x[0][0] = 0.5f * (x[0][1] + x[1][0]);
    x[NUM_ROW - 1][0] = 0.5f * (x[NUM_ROW - 1][1] + x[NUM_ROW - 2][0]);
    x[0][NUM_COL - 1] = 0.5f * (x[0][NUM_COL - 2] + x[1][NUM_COL - 1]);
    x[NUM_ROW - 1][NUM_COL - 1] = 0.5f * (x[NUM_ROW - 1][NUM_COL - 2] + x[NUM_ROW - 2][NUM_COL - 1]);
}

void forces(f3 u[NUM_ROW][NUM_COL], float timestep) {
    GUI_force* f;
    f3 impulse;
    f3 tmp;
    while (!queueEmpty()) {
        f = force_dequeue();
        impulse = f->force * 0.01;
        int x = (float)f->position.x;
        int y = (float)(f->position.y);
        u[y][x] = u[y][x] + impulse;
        free(f);
    }

    while (!ink_queueEmpty()) {
        f = ink_dequeue();
        impulse = f->force * timestep;
        for (int x = 1; x < NUM_COL - 1; x++) {
            for (int y = 1; y < NUM_ROW - 1; y++) {
                f3 pt = { (float)x, (float)y };
                float effect = (float)exp(-pt.dist(f->position) / f->impulseRadius);
                tmp = (effect * impulse);
                ink[y][x] = ink[y][x] + f3::abs(tmp);
            }
        }
        free(f);
    }
}

void advect(f3 f[NUM_ROW][NUM_COL], f3 x[NUM_ROW][NUM_COL], f3 f_prev[NUM_ROW][NUM_COL], float timestep, float gain, f3 borderConditions) {

    // inner elements
    for (int i = 1; i < NUM_COL - 1; i++) {
        for (int y = 1; y < NUM_ROW - 1; y++) {
            // follow the velocity field "back in time" 
            float x_old = (i - (x[y][i].x * timestep * NUM_COL)); if (x_old < 0.5) x_old = 0.5; if (x_old > NUM_COL - 1.5) x_old = NUM_COL - 1.5;
            float y_old = (y - (x[y][i].y * timestep * NUM_ROW)); if (y_old < 0.5) y_old = 0.5; if (y_old > NUM_ROW - 1.5) y_old = NUM_ROW - 1.5;


            int y_1 = (int)(y_old);
            int y_2 = y_1 + 1;
            int x_1 = (int)(x_old);
            int x_2 = x_1 + 1;


            // interpolate and update x
            float p_x = x_old - (int)x_old, p_y = y_old - (int)y_old;

            f[y][i] = (1 - p_y) * ((1 - p_x) * f_prev[y_1][x_1] + (p_x)*f_prev[y_1][x_2]) + (p_y) * ((1 - p_x) * f_prev[y_2][x_1] + (p_x)*f_prev[y_2][x_2]);
            f[y][i] = f[y][i] * gain;

        }
    }

    setBorders(f, borderConditions.x, borderConditions.y); //TODO: border condition
}

void diffuse(f3 x[NUM_ROW][NUM_COL], f3 x_prev[NUM_ROW][NUM_COL], float diffusion_rate, float dt, f3 borderConditions) {
    //compute inner elements

    float a = dt * diffusion_rate * GRID_SIZE;

    for (int iter = 0; iter < 20; iter++) {
        for (int i = 1; i < NUM_COL - 1; i++) {
            for (int j = 1; j < NUM_ROW - 1; j++) {
                x[j][i] = (x_prev[j][i] + a * (x[j - 1][i] + x[j + 1][i] + x[j][i - 1] + x[j][i + 1])) / (1 + 4 * a);
            }
        }
        setBorders(x, borderConditions.x, borderConditions.y);
    }
}

float divergence_u[NUM_ROW][NUM_COL];
void pressure(float p[NUM_ROW][NUM_COL], f3 u[NUM_ROW][NUM_COL], float timestep) {
    

    float h = 1.0 / NUM_COL;

    for (int y = 1; y < NUM_ROW - 1; y++) {
        for (int x = 1; x < NUM_COL - 1; x++) {
            divergence_u[y][x] = (u[y][x + 1].x - u[y][x - 1].x + u[y + 1][x].y - u[y - 1][x].y) * -0.5 * h; //dFx/dx
            p[y][x] = 0;
        }
    }
    setBorders(divergence_u, 0);
    setBorders(p, 0);

    for (int iter = 0; iter < 20; iter++) {
        for (int i = 1; i < NUM_COL - 1; i++) {
            for (int j = 1; j < NUM_ROW - 1; j++) {
                p[j][i] = (divergence_u[j][i] + p[j - 1][i] + p[j + 1][i] + p[j][i - 1] + p[j][i + 1]) * 0.25;
            }
        }
        setBorders(p, 0);
    }

    for (int y = 1; y < NUM_ROW - 1; y++) {
        for (int x = 1; x < NUM_COL - 1; x++) {
            u[y][x] = u[y][x] - 0.5 * f3(p[y][x + 1] - p[y][x - 1], p[y + 1][x] - p[y - 1][x]) * NUM_COL;
        }
    }
    setBorders(u, 1, 2);
}

void addExternalForce(f3 force, f3 position, float impulseRadius) {
    force_enqueue(force, position, impulseRadius);
}

void addInk(f3 force, f3 position, float impulseRadius) {
    ink_enqueue(force, position, impulseRadius);
}