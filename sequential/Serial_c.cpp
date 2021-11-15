#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "utilities.h"

#define min(x, y) (x > y) ? y : x
#define max(x, y) (x < y) ? y : x

#define NUM_ROW 5
#define NUM_COL 5
#define BOUNDARY_GAIN 0.75


enum initialGuess {
    ALL_ZEROS,
    USE_VAR
};

// update functions
void update_grid(f2*** u, f2*** u_temp, f2** p, f2** divergence_u, float timestep, float epsilon);
void advect(f2** u, f2** u1, float timestep);
void diffuse(f2** u, f2** u1);
void addForces(f2** u, float timestep);
void computePressure(f2** p, f2** u, f2** divergence_u);
void subtractPressureGradient(f2** u, f2** u1, f2** p);
void jacobi(f2** xNew, f2** guess, f2** b, float alpha, float rBeta, int numIterations, enum initialGuess use_guess);
void divergence(f2** div_f, f2** f);

// GUI event
void addExternalForce(f2 force, f2 position, float impulseRadius);

int main(int argc, char* argv[])
{
    float timestep;
    float epsilon;
    if (argc != 3 || argv[1] == NULL)
    {
        printf("ERROR:\tInput needs to be ./sequential <timestep> <epsilon>.\n");
        timestep = 0.1;
        epsilon = 0.1;
        //return 0;
    }
    else {
         timestep = (float)atof(argv[1]);
         epsilon = (float)atof(argv[2]);
    }

    

    if (timestep == 0)
    {
        printf("ERROR:\t<timestep> must be a positive integer larger than 0.\n");
        // return 0;
        timestep = 0.1;
    }

    if (epsilon == 0 || epsilon > 1)
    {
        printf("ERROR:\t<epsilon> must be a positive float in range (0, 1].\n");
        // return 0;
        epsilon = 0.001;
    }

    // assign pointers for u, u1, u2 2D grids
    f2** u, ** u_temp, ** p, ** divergence_u;

    // memory allocation of arrays
    u = (f2**)malloc(sizeof(f2*) * NUM_ROW);
    u_temp = (f2**)malloc(sizeof(f2*) * NUM_ROW);
    p = (f2**)malloc(sizeof(f2*) * NUM_ROW);
    divergence_u = (f2**)malloc(sizeof(f2*) * NUM_ROW);

    for (int j = 0; j < NUM_ROW; j++) {
        u[j] = (f2*)calloc(NUM_COL, sizeof(f2));
        u_temp[j] = (f2*)calloc(NUM_COL, sizeof(f2));
        p[j] = (f2*)calloc(NUM_COL, sizeof(f2));
        divergence_u[j] = (f2*)calloc(NUM_COL, sizeof(f2));
    }

    // simulate initial hit, update u by inrementing the timestep
    f2 f = { 20, 1 };
    f2 pos = { NUM_COL / 2, NUM_ROW / 2 };
    addExternalForce(f, pos, min(2, NUM_COL));

    //run through timesteps

    for (int q = 0; q < 10; q++) {
        update_grid(&u, &u_temp, p, divergence_u, timestep, epsilon);
        printf("done with iteration # %d\n", q);
        for (int i = 0; i < NUM_ROW; i++) {
            for (int j = 0; j < NUM_COL; j++) {
                printf("(%f,%f), ", u[i][j].x, u[i][j].y);
            }
            puts("");
        }
        puts("");
    }

    for (int j = 0; j < NUM_ROW; j++) {
        free(u[j]);
        free(u_temp[j]);
        free(p[j]);
        free(divergence_u[j]);
    }
    free(u);
    free(u_temp);
    free(p);
    free(divergence_u);

    return 0;
}

void update_grid(f2*** u, f2*** u_temp, f2** p, f2** divergence_u, float timestep, float epsilon) {

    // Apply the first 3 operators in Equation 12. 
    advect(*u_temp, *u, timestep);
    for (int i = 0; i < NUM_ROW; i++) {
        for (int j = 0; j < NUM_COL; j++) {
            printf("(%f,%f), ", (*u_temp)[i][j].x, (*u_temp)[i][j].y);
        }
        puts("");
    }
    puts("");
    //u1 will store the output of diffuse
    diffuse(*u, *u_temp);
    for (int i = 0; i < NUM_ROW; i++) {
        for (int j = 0; j < NUM_COL; j++) {
            printf("(%f,%f), ", (*u)[i][j].x, (*u)[i][j].y);
        }
        puts("");
    }
    puts("");
    addForces(*u, timestep);

    for (int i = 0; i < NUM_ROW; i++) {
        for (int j = 0; j < NUM_COL; j++) {
            printf("(%f,%f), ", (*u)[i][j].x, (*u)[i][j].y);
        }
        puts("");
    }
    puts("");
    divergence(divergence_u, *u);
    // Now apply the projection operator to the result. 
    computePressure(p, *u, divergence_u);

    subtractPressureGradient(*u_temp, *u, p);

    //swap u and u_temp
    f2** tmp = *u;
    *u = *u_temp;
    *u_temp = tmp;
}

//void advect(f2 coords : WPOS,   // grid coordinates     
//    out float4 xNew : COLOR,  // advected qty     
//    uniform float timestep,             
//    uniform    float rdx,        // 1 / grid scale     
//    uniform    samplerRECT u,    // input velocity     
//    uniform    samplerRECT x)    // qty to advect 
//{   
//    // follow the velocity field "back in time"     
//    f2 pos = coords - timestep * rdx * f2texRECT(u, coords); 
// // interpolate and write to the output fragment   
//    xNew = f4texRECTbilerp(x, pos); 
//} 
// u = input velocity; x = 'ink' being advected
void advect(f2** u, f2** u1, float timestep) {

    // inner elements
    for (int x = 1; x < NUM_COL - 1; x++) {
        for (int y = 1; y < NUM_ROW - 1; y++) {
            // follow the velocity field "back in time" 
            float x_old = min(NUM_COL - 1, max(0, x - u1[y][x].x * timestep));
            float y_old = min(NUM_ROW - 1, max(0, y - u1[y][x].y * timestep));
            int y_1 = min(NUM_ROW - 1, max(0, (int)(y_old - 0.5)));
            int y_2 = min(NUM_ROW - 1, max(0, (int)(y_old + 0.5)));
            int x_1 = min(NUM_COL - 1, max(0, (int)(x_old - 0.5)));
            int x_2 = min(NUM_COL - 1, max(0, (int)(x_old + 0.5)));
            // interpolate and update u
            u[y][x].x = (u1[y_1][x_1].x + u1[y_1][x_2].x + u1[y_2][x_1].x + u1[y_2][x_2].x) / (float)4.0;
            u[y][x].y = (u1[y_1][x_1].y + u1[y_1][x_2].y + u1[y_2][x_1].y + u1[y_2][x_2].y) / (float)4.0;
        }
    }

    // upper + lower border
    for (int x = 1; x < NUM_COL - 1; x++) {
        //upper border
        u[0][x].x = -u[1][x].x;
        u[0][x].y = -u[1][x].y;

        //lower border
        u[NUM_ROW - 1][x].x = -u[NUM_ROW - 2][x].x;
        u[NUM_ROW - 1][x].y = -u[NUM_ROW - 2][x].y;
    }

    //left + right border
    for (int y = 1; y < NUM_ROW - 1; y++) {
        //left border
        u[y][0].x = -u[y][1].x;
        u[y][0].y = -u[y][1].y;

        //right border
        u[y][NUM_COL - 1].x = -u[y][NUM_COL - 2].x;
        u[y][NUM_COL - 1].y = -u[y][NUM_COL - 2].y;
    }
}

//void jacobi(half2 coords : WPOS,   // grid coordinates     
//    out    half4 xNew : COLOR,  // result     
//    uniform    half alpha,             
//    uniform    half rBeta,      // reciprocal beta     
//    uniform samplerRECT x,   // x vector (Ax = b)     
//    uniform samplerRECT b)   // b vector (Ax = b) 
//{   // left, right, bottom, and top x samples    
//    half4 xL = h4texRECT(x, coords - half2(1, 0));   
//    half4 xR = h4texRECT(x, coords + half2(1, 0));   
//    half4 xB = h4texRECT(x, coords - half2(0, 1));   
//    half4 xT = h4texRECT(x, coords + half2(0, 1)); 
//    // b sample, from center     
//    half4 bC = h4texRECT(b, coords); 
//    // evaluate Jacobi iteration   
//    xNew = (xL + xR + xB + xT + alpha * bC) * rBeta; 
//} 
void jacobi(f2** xNew, f2** guess, f2** b, float alpha, float rBeta, int numIterations, enum initialGuess use_guess) {
    if (use_guess == ALL_ZEROS) {
        for (int x = 0; x < NUM_COL; x++) {
            for (int y = 0; y < NUM_ROW; y++) {
                xNew[y][x].x = alpha * b[y][x].x * rBeta;
                xNew[y][x].y = alpha * b[y][x].y * rBeta;
            }
        }
    }
    else {
        for (int x = 1; x < NUM_COL - 1; x++) {
            for (int y = 1; y < NUM_ROW - 1; y++) {
                xNew[y][x].x = (guess[y - 1][x].x + guess[y + 1][x].x + guess[y][x - 1].x + guess[y][x + 1].x + alpha * b[y][x].x) * rBeta;
                xNew[y][x].y = (guess[y - 1][x].y + guess[y + 1][x].y + guess[y][x - 1].y + guess[y][x + 1].y + alpha * b[y][x].y) * rBeta;
            }
        }
    }


    for (int iter = 1; iter < numIterations; iter++) {
        for (int x = 1; x < NUM_COL - 1; x++) {
            for (int y = 1; y < NUM_ROW - 1; y++) {
                xNew[y][x].x = (xNew[y - 1][x].x + xNew[y + 1][x].x + xNew[y][x - 1].x + xNew[y][x + 1].x + alpha * b[y][x].x) * rBeta;
                xNew[y][x].y = (xNew[y - 1][x].y + xNew[y + 1][x].y + xNew[y][x - 1].y + xNew[y][x + 1].y + alpha * b[y][x].y) * rBeta;
            }
        }
    }
}

void diffuse(f2** u, f2** u1) {
    //compute inner elements
    jacobi(u, u1, u1, 10, 10, 25, USE_VAR);

    // upper + lower border
    for (int x = 1; x < NUM_COL - 1; x++) {
        //upper border
        u[0][x].x = -u[1][x].x;
        u[0][x].y = -u[1][x].y;

        //lower border
        u[NUM_ROW - 1][x].x = -u[NUM_ROW - 2][x].x;
        u[NUM_ROW - 1][x].y = -u[NUM_ROW - 2][x].y;
    }

    //left + right border
    for (int y = 1; y < NUM_ROW - 1; y++) {
        //left border
        u[y][0].x = -u[y][1].x;
        u[y][0].y = -u[y][1].y;

        //right border
        u[y][NUM_COL - 1].x = -u[y][NUM_COL - 2].x;
        u[y][NUM_COL - 1].y = -u[y][NUM_COL - 2].y;
    }
}

void addForces(f2** u, float timestep) {
    GUI_force* f;
    f2 impulse;
    while (!queueEmpty()) {
        f = force_dequeue();
        impulse.x = f->force.x * timestep;
        impulse.y = f->force.y * timestep;

        for (int x = 0; x < NUM_COL; x++) {
            for (int y = 0; y < NUM_ROW; y++) {
                float effect = (float) exp(((x - f->position.x) * (x - f->position.x) + (y - f->position.y) * (y - f->position.y)) / f->impulseRadius);
                u[y][x].x = effect * impulse.x;
                u[y][x].y = effect * impulse.y;
            }
        }
        free(f);
    }
}

//void divergence(half2 coords : WPOS,   // grid coordinates     
//    out    half4 div : COLOR,  // divergence     
//    uniform half halfrdx,   // 0.5 / gridscale     
//    uniform samplerRECT w)  // vector field 
//{
//    half4 wL = h4texRECT(w, coords - half2(1, 0));
//    half4 wR = h4texRECT(w, coords + half2(1, 0));
//    half4 wB = h4texRECT(w, coords - half2(0, 1));
//    half4 wT = h4texRECT(w, coords + half2(0, 1));
//    div = halfrdx * ((wR.x - wL.x) + (wT.y - wB.y));
//}
void divergence(f2** div_f, f2** f) {
    // inner elements
    for (int x = 1; x < NUM_COL - 1; x++) {
        for (int y = 1; y < NUM_ROW - 1; y++) {
            div_f[y][x].x = (f[y][x + 1].x - f[y][x - 1].x) / (float)2.0;
            div_f[y][x].x = (f[y - 1][x].y - f[y + 1][x].y) / (float)2.0;
        }
    }

    // upper + lower border
    for (int x = 1; x < NUM_COL - 1; x++) {
        //upper border
        div_f[0][x].x = (f[1][x].x) / (float)2.0;
        div_f[0][x].y = (-f[1][x].y) / (float)2.0;

        //lower border
        div_f[NUM_ROW - 1][x].x = (-f[NUM_ROW - 2][x].x) / (float)2.0;
        div_f[NUM_ROW - 1][x].y = (f[NUM_ROW - 2][x].y) / (float)2.0;
    }

    //left + right border
    for (int y = 1; y < NUM_ROW - 1; y++) {
        //left border
        div_f[y][0].x = (f[y][1].x) / (float)2.0;
        div_f[y][0].y = (-f[y][1].y) / (float)2.0;

        //right border
        div_f[y][NUM_COL - 1].x = (-f[y][NUM_COL - 2].x) / (float)2.0;
        div_f[y][NUM_COL - 1].y = (f[y][NUM_COL - 2].y) / (float)2.0;
    }
}

void computePressure(f2** p, f2** u, f2** divergence_u) {
    jacobi(p, p, divergence_u, 10, 0.25, 25, ALL_ZEROS);

    // upper + lower border
    for (int x = 1; x < NUM_COL - 1; x++) {
        //upper border
        p[0][x].x = p[1][x].x;
        p[0][x].y = p[1][x].y;

        //lower border
        p[NUM_ROW - 1][x].x = p[NUM_ROW - 2][x].x;
        p[NUM_ROW - 1][x].y = p[NUM_ROW - 2][x].y;
    }

    //left + right border
    for (int y = 1; y < NUM_ROW - 1; y++) {
        //left border
        p[y][0].x = p[y][1].x;
        p[y][0].y = p[y][1].y;

        //right border
        p[y][NUM_COL - 1].x = p[y][NUM_COL - 2].x;
        p[y][NUM_COL - 1].y = p[y][NUM_COL - 2].y;
    }
}

//void gradient(half2 coords : WPOS,   // grid coordinates     
//    out half4 uNew : COLOR,  // new velocity     
//    uniform half halfrdx,    // 0.5 / gridscale     
//    uniform samplerRECT p,   // pressure     
//    uniform samplerRECT w)   // velocity 
//{   
//    half pL = h1texRECT(p, coords - half2(1, 0));   
//    half pR = h1texRECT(p, coords + half2(1, 0));   
//    half pB = h1texRECT(p, coords - half2(0, 1));   
//    half pT = h1texRECT(p, coords + half2(0, 1)); 
//    uNew = h4texRECT(w, coords);   
//    uNew.xy -= halfrdx * half2(pR - pL, pT - pB); 
//}
void subtractPressureGradient(f2** u, f2** u1, f2** p) {

}

/**
* TODO: Should add forces to a queue
*/
void addExternalForce(f2 force, f2 position, float impulseRadius) {
    force_enqueue(force, position, impulseRadius);
}