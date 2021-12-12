#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
// #include "freeglut.h"
#include "utilities.h"

// Defining static constant grid parameters
#define NUM_ROW 128
#define NUM_COL NUM_ROW
#define GRID_SIZE NUM_ROW *NUM_COL
#define BOUNDARY_GAIN 0.75

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Helper methods and classes for simulation.
std::ostream &operator<<(std::ostream &strm, const f3 &a)
{
    return strm << "f2(" << a.x << ", " << a.y << ", " << a.z << ")";
}

f3 &f3::operator=(const f3 &c1)
{
    x = c1.x;
    y = c1.y;
    z = c1.z;
    return *this;
} // f2

////////////////
// addition
// ////////////
f3 operator+(const f3 &c1, const f3 &c2) { return f3(c1.x + c2.x, c1.y + c2.y, c1.z + c2.z); } // f2
f3 operator+(const f3 &c1, const float c2) { return f3(c1.x + c2, c1.y + c2, c1.z + c2); }     // float
f3 operator+(const float c2, const f3 &c1) { return f3(c1.x + c2, c1.y + c2, c1.z + c2); }
f3 operator+(const f3 &c1, const double c2) { return f3(c1.x + (float)c2, c1.y + (float)c2, c1.z + (float)c2); } // double
f3 operator+(const double c2, const f3 &c1) { return f3(c1.x + (float)c2, c1.y + (float)c2, c1.z + (float)c2); }
f3 operator+(const f3 &c1, const int c2) { return f3(c1.x + c2, c1.y + c2, c1.z + c2); } // int
f3 operator+(const int c2, const f3 &c1) { return f3(c1.x + c2, c1.y + c2, c1.z + c2); }

////////////////
// subtraction
// ////////////
f3 operator-(const f3 &c1, const f3 &c2) { return f3(c1.x - c2.x, c1.y - c2.y, c1.z - c2.z); } // f2
f3 f3::operator-() const { return f3(-x, -y, -z); }
f3 operator-(const f3 &c1, const float c2) { return f3(c1.x - c2, c1.y - c2, c1.z - c2); } // float
f3 operator-(const float c2, const f3 &c1) { return f3(c1.x - c2, c1.y - c2, c1.z - c2); }
f3 operator-(const f3 &c1, const double c2) { return f3(c1.x - (float)c2, c1.y - (float)c2, c1.z - (float)c2); } // double
f3 operator-(const double c2, const f3 &c1) { return f3(c1.x - (float)c2, c1.y - (float)c2, c1.z - (float)c2); }
f3 operator-(const f3 &c1, const int c2) { return f3(c1.x - c2, c1.y - c2, c1.z - c2); } // int
f3 operator-(const int c2, const f3 &c1) { return f3(c1.x - c2, c1.y - c2, c1.z - c2); }

////////////////
// multiplication
// ////////////
f3 operator*(const f3 &c1, const f3 &c2) { return f3(c1.x * c2.x, c1.y * c2.y, c1.z * c2.z); } // f2
f3 operator*(const f3 &c1, const float c2) { return f3(c1.x * c2, c1.y * c2, c1.z * c2); }     // float
f3 operator*(const float c2, const f3 &c1) { return f3(c1.x * c2, c1.y * c2, c1.z * c2); }
f3 operator*(const f3 &c1, const double c2) { return f3(c1.x * (float)c2, c1.y * (float)c2, c1.z * (float)c2); } // double
f3 operator*(const double c2, const f3 &c1) { return f3(c1.x * (float)c2, c1.y * (float)c2, c1.z * (float)c2); }
f3 operator*(const f3 &c1, const int c2) { return f3(c1.x * c2, c1.y * c2, c1.z * c2); } // int
f3 operator*(const int c2, const f3 &c1) { return f3(c1.x * c2, c1.y * c2, c1.z * c2); }

////////////////
// division
// ////////////
f3 operator/(const f3 &c1, const f3 &c2) { return f3(c1.x / c2.x, c1.y / c2.y, c1.z / c2.z); } // f2
f3 operator/(const f3 &c1, const float c2) { return f3(c1.x / c2, c1.y / c2, c1.z / c2); }     // float
f3 operator/(const float c2, const f3 &c1) { return f3(c1.x / c2, c1.y / c2, c1.z / c2); }
f3 operator/(const f3 &c1, const double c2) { return f3(c1.x / (float)c2, c1.y / (float)c2, c1.z / (float)c2); } // double
f3 operator/(const double c2, const f3 &c1) { return f3(c1.x / (float)c2, c1.y / (float)c2, c1.z / (float)c2); }
f3 operator/(const f3 &c1, const int c2) { return f3(c1.x / c2, c1.y / c2, c1.z / c2); } // int
f3 operator/(const int c2, const f3 &c1) { return f3(c1.x / c2, c1.y / c2, c1.z / c2); }

queue forceQueue = {NULL, NULL, 0};
queue inkQueue = {NULL, NULL, 0};

extern void force_enqueue(f3 force, f3 position, float impulseRadius)
{
    node *newNode = (node *)malloc(sizeof(node));
    newNode->force = (GUI_force *)malloc(sizeof(GUI_force));
    newNode->force->force = force;
    newNode->force->position = position;
    newNode->force->impulseRadius = impulseRadius;
    newNode->next = NULL;
    newNode->previous = NULL;

    if (queueSize() == 0)
    {
        forceQueue.first = newNode;
        forceQueue.last = newNode;
    }
    else
    {
        forceQueue.last->next = newNode;
        newNode->previous = forceQueue.last;
        forceQueue.last = newNode;
    }
    forceQueue.size++;
}

extern GUI_force *force_dequeue()
{
    if (queueSize() == 0)
    {
        return NULL;
    }
    GUI_force *force;
    if (queueSize() == 1)
    {
        // get last element from queue
        force = forceQueue.last->force;

        // remove last element from queue, which is also the first element of queue
        free(forceQueue.last);
        forceQueue.first = NULL;
        forceQueue.last = NULL;
    }
    else
    {
        // get last element from queue
        force = forceQueue.last->force;

        // remove last element from queue
        node *toDelete = forceQueue.last;
        forceQueue.last = forceQueue.last->previous;
        forceQueue.last->next = NULL;
        free(toDelete);
    }
    forceQueue.size--;
    return force;
}

extern int queueEmpty()
{
    return queueSize() == 0;
}

extern int queueSize()
{
    return forceQueue.size;
}

extern void ink_enqueue(f3 force, f3 position, float impulseRadius)
{
    node *newNode = (node *)malloc(sizeof(node));
    newNode->force = (GUI_force *)malloc(sizeof(GUI_force));
    newNode->force->force = force;
    newNode->force->position = position;
    newNode->force->impulseRadius = impulseRadius;
    newNode->next = NULL;
    newNode->previous = NULL;

    if (ink_queueSize() == 0)
    {
        inkQueue.first = newNode;
        inkQueue.last = newNode;
    }
    else
    {
        inkQueue.last->next = newNode;
        newNode->previous = inkQueue.last;
        inkQueue.last = newNode;
    }
    inkQueue.size++;
}

extern GUI_force *ink_dequeue()
{
    if (ink_queueSize() == 0)
    {
        return NULL;
    }
    GUI_force *force;
    if (ink_queueSize() == 1)
    {
        // get last element from queue
        force = inkQueue.last->force;

        // remove last element from queue, which is also the first element of queue
        free(inkQueue.last);
        inkQueue.first = NULL;
        inkQueue.last = NULL;
    }
    else
    {
        // get last element from queue
        force = inkQueue.last->force;

        // remove last element from queue
        node *toDelete = inkQueue.last;
        inkQueue.last = inkQueue.last->previous;
        inkQueue.last->next = NULL;
        free(toDelete);
    }
    inkQueue.size--;
    return force;
}

extern int ink_queueEmpty()
{
    return ink_queueSize() == 0;
}

extern int ink_queueSize()
{
    return inkQueue.size;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// update functions
void forces(f3 u[NUM_ROW][NUM_COL], float timestep);
void advect(f3 f[NUM_ROW][NUM_COL], f3 x[NUM_ROW][NUM_COL], f3 f_prev[NUM_ROW][NUM_COL], float timestep, float gain, f3 borderConditions);
void diffuse(f3 x[NUM_ROW][NUM_COL], f3 x_prev[NUM_ROW][NUM_COL], float diffusion_rate, float dt, f3 borderConditions);
void pressure(float p[NUM_ROW][NUM_COL], f3 u[NUM_ROW][NUM_COL], float timestep);
void setBorders(f3 x[NUM_ROW][NUM_COL], int borderCondition);
void setBorders(float x[NUM_ROW][NUM_COL], int borderCondition);
void update_grid(float timestep, float epsilon);

// GUI event
void addExternalForce(f3 force, f3 position, float impulseRadius);
void addInk(f3 force, f3 position, float impulseRadius);

// assign pointers for u, u1, u2 2D grids
f3 U[NUM_ROW][NUM_COL], U_prev[NUM_ROW][NUM_COL];     // velocity fields
float P[NUM_ROW][NUM_COL], P_prev[NUM_ROW][NUM_COL];  // pressure fields
f3 ink[NUM_ROW][NUM_COL], ink_prev[NUM_ROW][NUM_COL]; // density fields

// Defining Global Simulation Parameters
volatile static float dt, diff, visc, vort, fluidSize, fluidAmount, forceSize, forceAmount;

// Clearing the grid.
void clearGrids()
{
    for (int x = 0; x < NUM_COL; x++)
    {
        for (int y = 0; y < NUM_ROW; y++)
        {
            P[y][x] = P_prev[y][x] = 0;
            U[y][x] = U_prev[y][x] = ink[y][x] = ink_prev[y][x] = f3(0, 0);
        }
    }
}

static void printFrameMatrices(f3 ink[NUM_ROW][NUM_COL], f3 U[NUM_ROW][NUM_COL], float P[NUM_ROW][NUM_COL], int N, int iteration)
{
    printf("Iteration: %d\n", iteration);
    printf("Density Matrix:\n");
    for (int i = 0; i <= N + 1; i++)
    {
        for (int j = 0; j <= N + 1; j++)
        {
            std::cout << ink[i][j] << ", ";
        }
        printf("\n");
    }
    printf("\n");
    printf("Velocity U Matrix:\n");
    for (int i = 0; i <= N + 1; i++)
    {
        for (int j = 0; j <= N + 1; j++)
        {
            std::cout << U[i][j] << ", ";
        }
        printf("\n");
    }
    printf("\n");
    printf("Pressure P Matrix:\n");
    for (int i = 0; i <= N + 1; i++)
    {
        for (int j = 0; j <= N + 1; j++)
        {
            std::cout << P[i][j] << ", ";
        }
        printf("\n");
    }
    printf("\n");
}

/*
  ----------------------------------------------------------------------
   main --- main routine
  ----------------------------------------------------------------------
*/

int main(int argc, char **argv)
{
    // declare variables
    dt = 0.1f;
    diff = 0.1f;
    visc = 0.1f;
    float force = 100.0f, source = 100.0f;
    int N = 128, num_iterations = 10, display_output = 0;

    // load parameters from command line
    if (argc != 3)
    {
        fprintf(stderr, "usage : %s <number of iterations> <display output>\n", argv[0]);
        fprintf(stderr, "where:\n");
        fprintf(stderr, "\t number_of_iterations     : number of iterations to run the simulation\n");
        fprintf(stderr, "\t display_output           : Whether to display the opdated grids after every iteration\n");
        exit(1);
    }

    num_iterations = atoi(argv[1]);
    display_output = atoi(argv[2]);

    // validate parameters
    if (num_iterations <= 0)
    {
        fprintf(stderr, "number_of_iterations must be an integer larger than 0\n");
        exit(1);
    }
    if (display_output != 0 && display_output != 1)
    {
        fprintf(stderr, "display_output must be 1 if TRUE and 0 if FALSE\n");
        exit(1);
    }

    // Setting added vars
    dt = 0.1;
    visc = 0.01;
    vort = 0.5;
    fluidSize = 3;
    fluidAmount = 0.1;
    forceSize = 2;
    forceAmount = 7;

    // simulate initial hit, update u by incrementing the timestep
    clearGrids();
    f3 pos = {NUM_COL / 2, NUM_ROW / 2};
    f3 f = {fluidAmount, fluidAmount};
    addInk(f, pos, fluidSize);
    f = {forceAmount, forceAmount};
    addExternalForce(f, pos, forceSize);

    if (display_output == 1)
    {
        printFrameMatrices(ink, U, P, NUM_ROW, 0);
    }

    // update loop
    clock_t start = clock();
    for (int i = 1; i <= num_iterations; i++)
    {
        update_grid(dt, visc);
        if (display_output == 1)
        {
            printFrameMatrices(ink, U, P, NUM_ROW, i);
        }
    }
    clock_t stop = clock();

    double elapsed = (double)(stop - start) * 1000.0 / CLOCKS_PER_SEC;
    printf("Runtime in ms: %f \n", elapsed);

    exit(0);
}

void update_grid(float timestep, float epsilon)
{

    forces(U, timestep);

    // velocity update
    diffuse(U_prev, U, epsilon / 2, timestep, f3(1, 2));
    pressure(P, U_prev, timestep);
    advect(U, U_prev, U_prev, timestep, 1, f3(1, 2));
    pressure(P, U, timestep);

    // density step
    diffuse(ink_prev, ink, 0.0001, timestep, f3(3, 3)); // TODO: have dynamic memory allocation
    advect(ink, U, ink_prev, timestep, 1, f3(3, 3));
}

void setBorders(f3 x[NUM_ROW][NUM_COL], int borderCondition)
{
    // upper + lower border
    for (int ixd = 0; ixd < NUM_COL; ixd++)
    {
        switch (borderCondition)
        {
        case 1:
            x[0][ixd] = -x[1][ixd];                     // upper border
            x[NUM_ROW - 1][ixd] = -x[NUM_ROW - 2][ixd]; // lower border
            x[ixd][0] = x[ixd][1];                      // left border
            x[ixd][NUM_COL - 1] = x[ixd][NUM_COL - 2];  // right border
            break;
        case 2:
            x[0][ixd] = x[1][ixd];                      // upper border
            x[NUM_ROW - 1][ixd] = x[NUM_ROW - 2][ixd];  // lower border
            x[ixd][0] = -x[ixd][1];                     // left border
            x[ixd][NUM_COL - 1] = -x[ixd][NUM_COL - 2]; // right border
            break;
        default:
            x[0][ixd] = x[1][ixd];                     // upper border
            x[NUM_ROW - 1][ixd] = x[NUM_ROW - 2][ixd]; // lower border
            x[ixd][0] = x[ixd][1];                     // left border
            x[ixd][NUM_COL - 1] = x[ixd][NUM_COL - 2]; // right border
            break;
        }
    }

    x[0][0] = 0.5f * (x[0][1] + x[1][0]);
    x[NUM_ROW - 1][0] = 0.5f * (x[NUM_ROW - 1][1] + x[NUM_ROW - 2][0]);
    x[0][NUM_COL - 1] = 0.5f * (x[0][NUM_COL - 2] + x[1][NUM_COL - 1]);
    x[NUM_ROW - 1][NUM_COL - 1] = 0.5f * (x[NUM_ROW - 1][NUM_COL - 2] + x[NUM_ROW - 2][NUM_COL - 1]);
}
void setBorders(f3 x[NUM_ROW][NUM_COL], int borderConditionX, int borderConditionY)
{
    // upper + lower border
    for (int ixd = 0; ixd < NUM_COL; ixd++)
    {
        switch (borderConditionX)
        {
        case 1:
            x[0][ixd].x = -x[1][ixd].x;                     // upper border
            x[NUM_ROW - 1][ixd].x = -x[NUM_ROW - 2][ixd].x; // lower border
            x[ixd][0].x = x[ixd][1].x;                      // left border
            x[ixd][NUM_COL - 1].x = x[ixd][NUM_COL - 2].x;  // right border
            break;
        case 2:
            x[0][ixd].x = x[1][ixd].x;                      // upper border
            x[NUM_ROW - 1][ixd].x = x[NUM_ROW - 2][ixd].x;  // lower border
            x[ixd][0].x = -x[ixd][1].x;                     // left border
            x[ixd][NUM_COL - 1].x = -x[ixd][NUM_COL - 2].x; // right border
            break;
        default:
            x[0][ixd].x = x[1][ixd].x;                     // upper border
            x[NUM_ROW - 1][ixd].x = x[NUM_ROW - 2][ixd].x; // lower border
            x[ixd][0].x = x[ixd][1].x;                     // left border
            x[ixd][NUM_COL - 1].x = x[ixd][NUM_COL - 2].x; // right border
            break;
        }
        switch (borderConditionY)
        {
        case 1:
            x[0][ixd].y = -x[1][ixd].y;                     // upper border
            x[NUM_ROW - 1][ixd].y = -x[NUM_ROW - 2][ixd].y; // lower border
            x[ixd][0].y = x[ixd][1].y;                      // left border
            x[ixd][NUM_COL - 1].y = x[ixd][NUM_COL - 2].y;  // right border
            break;
        case 2:
            x[0][ixd].y = x[1][ixd].y;                      // upper border
            x[NUM_ROW - 1][ixd].y = x[NUM_ROW - 2][ixd].y;  // lower border
            x[ixd][0].y = -x[ixd][1].y;                     // left border
            x[ixd][NUM_COL - 1].y = -x[ixd][NUM_COL - 2].y; // right border
            break;
        default:
            x[0][ixd].y = x[1][ixd].y;                     // upper border
            x[NUM_ROW - 1][ixd].y = x[NUM_ROW - 2][ixd].y; // lower border
            x[ixd][0].y = x[ixd][1].y;                     // left border
            x[ixd][NUM_COL - 1].y = x[ixd][NUM_COL - 2].y; // right border
            break;
        }
    }

    x[0][0] = 0.5f * (x[0][1] + x[1][0]);
    x[NUM_ROW - 1][0] = 0.5f * (x[NUM_ROW - 1][1] + x[NUM_ROW - 2][0]);
    x[0][NUM_COL - 1] = 0.5f * (x[0][NUM_COL - 2] + x[1][NUM_COL - 1]);
    x[NUM_ROW - 1][NUM_COL - 1] = 0.5f * (x[NUM_ROW - 1][NUM_COL - 2] + x[NUM_ROW - 2][NUM_COL - 1]);
}
void setBorders(float x[NUM_ROW][NUM_COL], int borderCondition)
{
    // upper + lower border
    for (int ixd = 0; ixd < NUM_COL; ixd++)
    {
        switch (borderCondition)
        {
        case 1:
            x[0][ixd] = -x[1][ixd];                     // upper border
            x[NUM_ROW - 1][ixd] = -x[NUM_ROW - 2][ixd]; // lower border
            x[ixd][0] = x[ixd][1];                      // left border
            x[ixd][NUM_COL - 1] = x[ixd][NUM_COL - 2];  // right border
            break;
        case 2:
            x[0][ixd] = x[1][ixd];                      // upper border
            x[NUM_ROW - 1][ixd] = x[NUM_ROW - 2][ixd];  // lower border
            x[ixd][0] = -x[ixd][1];                     // left border
            x[ixd][NUM_COL - 1] = -x[ixd][NUM_COL - 2]; // right border
            break;
        default:
            x[0][ixd] = x[1][ixd];                     // upper border
            x[NUM_ROW - 1][ixd] = x[NUM_ROW - 2][ixd]; // lower border
            x[ixd][0] = x[ixd][1];                     // left border
            x[ixd][NUM_COL - 1] = x[ixd][NUM_COL - 2]; // right border
            break;
        }
    }

    x[0][0] = 0.5f * (x[0][1] + x[1][0]);
    x[NUM_ROW - 1][0] = 0.5f * (x[NUM_ROW - 1][1] + x[NUM_ROW - 2][0]);
    x[0][NUM_COL - 1] = 0.5f * (x[0][NUM_COL - 2] + x[1][NUM_COL - 1]);
    x[NUM_ROW - 1][NUM_COL - 1] = 0.5f * (x[NUM_ROW - 1][NUM_COL - 2] + x[NUM_ROW - 2][NUM_COL - 1]);
}

void forces(f3 u[NUM_ROW][NUM_COL], float timestep)
{
    GUI_force *f;
    f3 impulse;
    f3 tmp;
    while (!queueEmpty())
    {
        f = force_dequeue();
        impulse = f->force * 0.01;
        int x = (float)f->position.x;
        int y = (float)(f->position.y);
        // u[y][x] = u[y][x] + impulse;
        for (int x = 1; x < NUM_COL - 1; x++)
        {
            for (int y = 1; y < NUM_ROW - 1; y++)
            {
                f3 pt = {(float)x, (float)y};

                float effect = (float)exp(-pt.dist(f->position) / f->impulseRadius);
                tmp = (effect * impulse);
                u[y][x] = u[y][x] + tmp;
            }
        }
        free(f);
    }

    int inkUpdated = 0;
    while (!ink_queueEmpty())
    {
        inkUpdated = 1;
        f = ink_dequeue();
        impulse = f->force * timestep;
        for (int x = 1; x < NUM_COL - 1; x++)
        {
            for (int y = 1; y < NUM_ROW - 1; y++)
            {
                f3 pt = {(float)x, (float)y};
                float effect = (float)exp(-pt.dist(f->position) / f->impulseRadius);
                tmp = (effect * impulse);
                ink[y][x] = ink[y][x] + f3::abs(tmp);
            }
        }
        free(f);
    }
    if (inkUpdated)
    {
        // TODO: SWAP inks?
    }
}
void advect(f3 f[NUM_ROW][NUM_COL], f3 x[NUM_ROW][NUM_COL], f3 f_prev[NUM_ROW][NUM_COL], float timestep, float gain, f3 borderConditions)
{

    // inner elements
    for (int i = 1; i < NUM_COL - 1; i++)
    {
        for (int y = 1; y < NUM_ROW - 1; y++)
        {
            // follow the velocity field "back in time"
            float x_old = (i - (x[y][i].x * timestep * NUM_COL));
            if (x_old < 0.5)
                x_old = 0.5;
            if (x_old > NUM_COL - 1.5)
                x_old = NUM_COL - 1.5;
            float y_old = (y - (x[y][i].y * timestep * NUM_ROW));
            if (y_old < 0.5)
                y_old = 0.5;
            if (y_old > NUM_ROW - 1.5)
                y_old = NUM_ROW - 1.5;

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

    setBorders(f, borderConditions.x, borderConditions.y); // TODO: border condition
}
void diffuse(f3 x[NUM_ROW][NUM_COL], f3 x_prev[NUM_ROW][NUM_COL], float diffusion_rate, float dt, f3 borderConditions)
{
    // compute inner elements

    float a = dt * diffusion_rate * GRID_SIZE;

    for (int iter = 0; iter < 20; iter++)
    {
        for (int i = 1; i < NUM_COL - 1; i++)
        {
            for (int j = 1; j < NUM_ROW - 1; j++)
            {
                x[j][i] = (x_prev[j][i] + a * (x[j - 1][i] + x[j + 1][i] + x[j][i - 1] + x[j][i + 1])) / (1 + 4 * a);
            }
        }
        setBorders(x, borderConditions.x, borderConditions.y); // TODO: boundary condition
    }
}

float divergence_u[NUM_ROW][NUM_COL];
void pressure(float p[NUM_ROW][NUM_COL], f3 u[NUM_ROW][NUM_COL], float timestep)
{

    float h = 1.0 / NUM_COL;

    for (int y = 1; y < NUM_ROW - 1; y++)
    {
        for (int x = 1; x < NUM_COL - 1; x++)
        {
            divergence_u[y][x] = (u[y][x + 1].x - u[y][x - 1].x + u[y + 1][x].y - u[y - 1][x].y) * -0.5 * h; // dFx/dx
            p[y][x] = 0;
        }
    }
    setBorders(divergence_u, 0);
    setBorders(p, 0);

    for (int iter = 0; iter < 20; iter++)
    {
        for (int i = 1; i < NUM_COL - 1; i++)
        {
            for (int j = 1; j < NUM_ROW - 1; j++)
            {
                p[j][i] = (divergence_u[j][i] + p[j - 1][i] + p[j + 1][i] + p[j][i - 1] + p[j][i + 1]) * 0.25;
            }
        }
        setBorders(p, 0);
    }

    for (int y = 1; y < NUM_ROW - 1; y++)
    {
        for (int x = 1; x < NUM_COL - 1; x++)
        {
            u[y][x] = u[y][x] - 0.5 * f3(p[y][x + 1] - p[y][x - 1], p[y + 1][x] - p[y - 1][x]) * NUM_COL;
        }
    }
    setBorders(u, 1, 2);
}

void curl(float curl_f[NUM_ROW][NUM_COL], f3 f[NUM_ROW][NUM_COL])
{
    // inner elements
    for (int x = 1; x < NUM_COL - 1; x++)
    {
        for (int y = 1; y < NUM_ROW - 1; y++)
        {
            curl_f[y][x] = (f[y][x + 1].y - f[y][x - 1].y) / (float)(2.0) - (f[y + 1][x].x - f[y - 1][x].x) / (float)(2.0);
        }
    }

    // TODO: fix
    //  upper + lower border
    for (int x = 0; x < NUM_COL; x++)
    {
        // upper border
        curl_f[0][x] = 0; // (f[1][x].x) / (float)2.0;

        // lower border
        curl_f[NUM_ROW - 1][x] = 0; //(-f[NUM_ROW - 2][x].x) / (float)2.0;
    }

    // left + right border
    for (int y = 0; y < NUM_ROW; y++)
    {
        // left border
        curl_f[y][0] = 0; //(f[y][1].x) / (float)2.0;

        // right border
        curl_f[y][NUM_COL - 1] = 0; //(-f[y][NUM_COL - 2].x) / (float)2.0;
    }
}

void gradient(f3 grad_f[NUM_ROW][NUM_COL], float f[NUM_ROW][NUM_COL])
{
    // inner elements
    for (int x = 1; x < NUM_COL - 1; x++)
    {
        for (int y = 1; y < NUM_ROW - 1; y++)
        {
            grad_f[y][x].x = (f[y][x + 1] - f[y][x - 1]) / (float)(2.0); // dFx/dx
            grad_f[y][x].y = (f[y + 1][x] - f[y - 1][x]) / (float)(2.0); // dFy/dy
        }
    }

    // TODO:fixe
    //  upper + lower border
    for (int x = 0; x < NUM_COL; x++)
    {
        // upper border
        grad_f[0][x].x = 0; //-grad_f[1][x].x;
        grad_f[0][x].y = 0; // grad_f[1][x].y;

        // lower border
        grad_f[NUM_ROW - 1][x].x = 0; //-grad_f[NUM_ROW - 2][x].x;
        grad_f[NUM_ROW - 1][x].y = 0; // grad_f[NUM_ROW - 2][x].y;
    }

    // left + right border
    for (int y = 0; y < NUM_ROW; y++)
    {
        // left border
        grad_f[y][0].x = 0; // grad_f[y][1].x;
        grad_f[y][0].y = 0; //-grad_f[y][1].y;

        // right border
        grad_f[y][NUM_COL - 1].x = 0; // grad_f[y][NUM_COL - 2].x;
        grad_f[y][NUM_COL - 1].y = 0; //-grad_f[y][NUM_COL - 2].y;
    }

    grad_f[0][0] = 0.5f * (grad_f[0][1] + grad_f[1][0]);
    grad_f[NUM_ROW - 1][0] = 0.5f * (grad_f[NUM_ROW - 1][1] + grad_f[NUM_ROW - 2][0]);
    grad_f[0][NUM_COL - 1] = 0.5f * (grad_f[0][NUM_COL - 2] + grad_f[1][NUM_COL - 1]);
    grad_f[NUM_ROW - 1][NUM_COL - 1] = 0.5f * (grad_f[NUM_ROW - 1][NUM_COL - 2] + grad_f[NUM_ROW - 2][NUM_COL - 1]);
}

void addToMatrix(f3 a[NUM_ROW][NUM_COL], f3 b[NUM_ROW][NUM_COL])
{
    for (int x = 0; x < NUM_COL; x++)
    {
        for (int y = 0; y < NUM_ROW; y++)
        {
            a[y][x] = a[y][x] + b[y][x];
        }
    }
}

void normalizeMatrix(f3 f[NUM_ROW][NUM_COL], float scalingFactor)
{
    for (int x = 0; x < NUM_COL; x++)
    {
        for (int y = 0; y < NUM_ROW; y++)
        {
            if (f[y][x].x != 0)
                f[y][x].x = f[y][x].x / abs(f[y][x].x) * scalingFactor;
            if (f[y][x].y != 0)
                f[y][x].y = f[y][x].y / abs(f[y][x].y) * scalingFactor;
            f[y][x].cap(0.1);
        }
    }
}

void matrixMultiply(f3 f[NUM_ROW][NUM_COL], float f2[NUM_ROW][NUM_COL])
{
    for (int x = 0; x < NUM_COL; x++)
    {
        for (int y = 0; y < NUM_ROW; y++)
        {
            if (abs(f2[y][x]) < 1e-5)
            {
                f[y][x] = {0, 0};
            }
        }
    }
}

void addExternalForce(f3 force, f3 position, float impulseRadius)
{
    force_enqueue(force, position, impulseRadius);
}

void addInk(f3 force, f3 position, float impulseRadius)
{
    ink_enqueue(force, position, impulseRadius);
}
