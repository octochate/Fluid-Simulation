#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "freeglut.h"
#include "utilities.h"




#define NUM_ROW 64
#define NUM_COL NUM_ROW
#define GRID_SIZE NUM_ROW * NUM_COL
#define BOUNDARY_GAIN 0.75

static void get_from_UI();

// update functions
void update_grid(f2 u[NUM_ROW][NUM_COL], f2 u_temp[NUM_ROW][NUM_COL], f2 p[NUM_ROW][NUM_COL], f2 divergence_u[NUM_ROW][NUM_COL], float timestep, float epsilon);
void advect(f2 u[NUM_ROW][NUM_COL], f2 u1[NUM_ROW][NUM_COL], f2 X[NUM_ROW][NUM_COL], float timestep);
void diffuse(f2 u[NUM_ROW][NUM_COL], f2 u1[NUM_ROW][NUM_COL], float diffusion_rate, float timestep);
void addForces(f2 u[NUM_ROW][NUM_COL], float timestep);
void computePressure(f2 p[NUM_ROW][NUM_COL], f2 u[NUM_ROW][NUM_COL], f2 divergence_u[NUM_ROW][NUM_COL], float timestep);
void subtractPressureGradient(f2 u[NUM_ROW][NUM_COL], f2 u1[NUM_ROW][NUM_COL], f2 p[NUM_ROW][NUM_COL]);
void divergence(f2 div_f[NUM_ROW][NUM_COL], f2 f[NUM_ROW][NUM_COL]);

// GUI event
void addExternalForce(f2 force, f2 position, float impulseRadius);

// assign pointers for u, u1, u2 2D grids
f2 U[2][NUM_ROW][NUM_COL], P[NUM_ROW][NUM_COL], Divergence_u[NUM_ROW][NUM_COL], ink[2][NUM_ROW][NUM_COL];

volatile float timestep = 0.01;
volatile float epsilon = 0.001;
volatile int q = 0;
volatile int q_ink = 0;

void clearGrids() {
    for (int x = 0; x < NUM_COL; x++) {
        for (int y = 0; y < NUM_ROW; y++) {
            P[y][x] = Divergence_u[y][x] = U[0][y][x] = U[1][y][x] = ink[0][y][x] = ink[1][y][x] = f2(0,0);
        }
    }
}

void RenderString(const char* string, float x, float y)
{
	glColor3f(1.0, 1.0, 1.0);
	glRasterPos2f(x, y);
}

static int win_id;
static int win_x, win_y;
static float dt, diff, visc;

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
    update_grid(U[q % 1], U[(q + 1) % 1], P, Divergence_u, timestep, epsilon);
    q = (q + 1) % 100;


	glutSetWindow(win_id);
	glutPostRedisplay();
}

static void display_func(void)
{
	pre_display();

	draw_velocity();

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

    /*if (i<1 || i>N || j<1 || j>N) return;

    if (mouse_down[0]) {
        u[IX(i, j)] = force * (mx - omx);
        v[IX(i, j)] = force * (omy - my);
    }

    if (mouse_down[2]) {
        d[IX(i, j)] = source;
    }*/

    f2 f = { 0.02f * (mx - omx) / timestep, 0.02f * (omy - my) / timestep };
    f2 pos =  {i, j};//{ NUM_COL / 2, NUM_ROW / 2 }; //
    addExternalForce(f, pos, 0.5);
    puts("mouses");

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

    // simulate initial hit, update u by inrementing the timestep
    

    clearGrids();

	win_x = 512;
	win_y = 512;
	open_glut_window();

	glutMainLoop();

	exit(0);
}


void update_grid(f2 u[NUM_ROW][NUM_COL], f2 u_temp[NUM_ROW][NUM_COL], f2 p[NUM_ROW][NUM_COL], f2 divergence_u[NUM_ROW][NUM_COL], float timestep, float epsilon) {

    // Apply the first 3 operators in Equation 12. 
    advect(u_temp, u, u, timestep);
    advect(ink[(q_ink + 1) % 1], u, ink[q_ink % 1], timestep);
    q_ink = (q_ink + 1) % 1;

    //u1 will store the output of diffuse
    diffuse(u, u_temp, epsilon, timestep);

    addForces(u, timestep);


    divergence(divergence_u, u);

 
    computePressure(p, u, divergence_u, timestep);

    subtractPressureGradient(u_temp, u, p);
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
void advect(f2 u[NUM_ROW][NUM_COL], f2 u1[NUM_ROW][NUM_COL], f2 X[NUM_ROW][NUM_COL], float timestep) {

    // inner elements
    for (int x = 0; x < NUM_COL; x++) {
        for (int y = 0; y < NUM_ROW; y++) {
            // follow the velocity field "back in time" 
            float x_old = (int)(x - (u1[y][x].x * timestep) + NUM_COL) % NUM_COL + 0.5;
            float y_old = (int)(y - (u1[y][x].y * timestep) + NUM_ROW) % NUM_ROW + 0.5;
            int y_1 = (int)(y_old - 0.5 + NUM_ROW) % NUM_ROW;
            int y_2 = (int)(y_old + 0.5) % NUM_ROW;
            int x_1 = (int)(x_old - 0.5 + NUM_COL) % NUM_COL;
            int x_2 = (int)(x_old + 0.5) % NUM_COL;
            // interpolate and update u
            u[y][x] = (X[y_1][x_1] + X[y_1][x_2] + X[y_2][x_1] + X[y_2][x_2]) / 4.0;
        }
    }

    // upper + lower border
    for (int x = 0; x < NUM_COL; x++) {
        //upper border
        u[0][x] =  -u[1][x];

        //lower border
        u[NUM_ROW - 1][x] = -u[NUM_ROW - 2][x];
    }

    //left + right border
    for (int y = 0; y < NUM_ROW; y++) {
        //left border
        u[y][0] = -u[y][1];

        //right border
        u[y][NUM_COL - 1] = -u[y][NUM_COL - 2];
    }

    // corners

}


void jacobi_diffuse(f2 xNew[NUM_ROW][NUM_COL], f2 guess[NUM_ROW][NUM_COL], f2 b[NUM_ROW][NUM_COL], int numIterations, float dt) {

    for (int x = 1; x < NUM_COL - 1; x++) {
        for (int y = 1; y < NUM_ROW - 1; y++) {
            f2 alpha = (guess[y][x] * guess[y][x] / dt);
            xNew[y][x] = (guess[y - 1][x] + guess[y + 1][x] + guess[y][x - 1] + guess[y][x + 1] + alpha * b[y][x]) * 1 / (4 + alpha);
        }
    }
 


    for (int iter = 1; iter < numIterations; iter++) {
        for (int x = 1; x < NUM_COL - 1; x++) {
            for (int y = 1; y < NUM_ROW - 1; y++) {
                f2 alpha = (guess[y][x] * guess[y][x] / dt);
                xNew[y][x] = (xNew[y - 1][x] + xNew[y + 1][x] + xNew[y][x - 1] + xNew[y][x + 1] + alpha * b[y][x]) * 1 / (4 + alpha);
            }
        }
    }
}

void diffuse(f2 u[NUM_ROW][NUM_COL], f2 u1[NUM_ROW][NUM_COL], float diffusion_rate, float timestep) {

    //compute inner elements
    jacobi_diffuse(u, u1, u1, 50, timestep * diffusion_rate);

    // upper + lower border
    for (int x = 0; x < NUM_COL; x++) {
        //upper border
        u[0][x] = -u[1][x];

        //lower border
        u[NUM_ROW - 1][x] = -u[NUM_ROW - 2][x];
    }

    //left + right border
    for (int y = 0; y < NUM_ROW; y++) {
        //left border
        u[y][0] = -u[y][1];

        //right border
        u[y][NUM_COL - 1] = -u[y][NUM_COL - 2];
    }
}

void addForces(f2 u[NUM_ROW][NUM_COL], float timestep) {
    GUI_force* f;
    f2 impulse;
    int inkUpdated = 0;
    while (!queueEmpty()) {
        inkUpdated = 1;
        f = force_dequeue();
        impulse = f->force * timestep;

        for (int x = 1; x < NUM_COL - 1; x++) {
            for (int y = 1; y < NUM_ROW - 1; y++) {
                float effect = (float)exp(-((x - f->position.x) * (x - f->position.x) + (y - f->position.y) * (y - f->position.y)) / f->impulseRadius);
                u[y][x] = u[y][x] + (effect * impulse);
                ink[(q_ink) % 1][y][x] = ink[(q_ink) % 1][y][x] + (effect * impulse);
                ink[(q_ink) % 1][y][x].cap(-1000, 1000);
            }
        }
        free(f);
    }

    if (inkUpdated) {
        q_ink = (q_ink + 1) % 1;
    }
}

void divergence(f2 div_f[NUM_ROW][NUM_COL], f2 f[NUM_ROW][NUM_COL]) {
    // inner elements
    for (int x = 1; x < NUM_COL; x++) {
        for (int y = 1; y < NUM_ROW; y++) {
            div_f[y][x].x = (f[y][x + 1].x - f[y][x - 1].x) / (float)2.0;
            div_f[y][x].x = (f[y - 1][x].y - f[y + 1][x].y) / (float)2.0;
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

void jacobi_pressure(f2 xNew[NUM_ROW][NUM_COL], f2 X[NUM_ROW][NUM_COL], f2 b[NUM_ROW][NUM_COL], int numIterations, float dt) {
    for (int x = 0; x < NUM_COL; x++) {
        for (int y = 0; y < NUM_ROW; y++) {
            xNew[y][x] = f2(0, 0);
        }
    }


    for (int iter = 1; iter < numIterations; iter++) {
        for (int x = 1; x < NUM_COL - 1; x++) {
            for (int y = 1; y < NUM_ROW - 1; y++) {
                f2 alpha = -X[y][x] * X[y][x];
                xNew[y][x] = (xNew[y - 1][x] + xNew[y + 1][x] + xNew[y][x - 1] + xNew[y][x + 1] + alpha * b[y][x]) * 0.25;
            }
        }
    }
}


void computePressure(f2 p[NUM_ROW][NUM_COL], f2 u[NUM_ROW][NUM_COL], f2 divergence_u[NUM_ROW][NUM_COL], float timestep) {
    jacobi_pressure(p, u, divergence_u, 40, timestep);

    // upper + lower border
    for (int x = 0; x < NUM_COL; x++) {
        //upper border
        p[0][x] = p[1][x];

        //lower border
        p[NUM_ROW - 1][x] = p[NUM_ROW - 2][x];
    }

    //left + right border
    for (int y = 0; y < NUM_ROW; y++) {
        //left border
        p[y][0] = p[y][1];

        //right border
        p[y][NUM_COL - 1] = p[y][NUM_COL - 2];
    }
}


void subtractPressureGradient(f2 u[NUM_ROW][NUM_COL], f2 u1[NUM_ROW][NUM_COL], f2 p[NUM_ROW][NUM_COL]) {
    for (int y = 1; y < NUM_ROW - 1; y++) {
        for (int x = 1; x < NUM_COL - 1; x++) {
            u[y][x] = u1[y][x] - f2(p[y][x + 1].x - p[y][x - 1].x, p[y - 1][x].y - p[y + 1][x].y);
        }
    }
}

/**
* TODO: Should add forces to a queue
*/
void addExternalForce(f2 force, f2 position, float impulseRadius) {
    force_enqueue(force, position, impulseRadius);
}