#include <stdlib.h>
#include <stdio.h>

#define IX(i,j) ((i)+(N+2)*(j))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}

void add_source(int N, float *x, float *s, float dt)
{
    int i, size = (N + 2) * (N + 2);
    for (i = 0; i < size; i++)
        x[i] += dt * s[i];
}

void set_bnd(int N, int b, float *x)
{
    int i;

    for (i = 1; i <= N; i++)
    {
        x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
        x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
    }
    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N + 1)] = 0.5f * (x[IX(1, N + 1)] + x[IX(0, N)]);
    x[IX(N + 1, 0)] = 0.5f * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
    x[IX(N + 1, N + 1)] = 0.5f * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}

void lin_solve(int N, int b, float *x, float *x0, float a, float c)
{
    int i, j, k;

    for (k = 0; k < 20; k++)
    {
        FOR_EACH_CELL
        x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] + x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
        END_FOR
        set_bnd(N, b, x);
    }
}

void diffuse(int N, int b, float *x, float *x0, float diff, float dt)
{
    float a = dt * diff * N * N;
    lin_solve(N, b, x, x0, a, 1 + 4 * a);
}

void advect ( int N, int b, float * d, float * d0, float * u, float * v, float dt )
{
	int i, j, i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1, dt0;

	dt0 = dt*N;
	FOR_EACH_CELL
		x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)];
		if (x<0.5f) x=0.5f; if (x>N+0.5f) x=N+0.5f; i0=(int)x; i1=i0+1;
		if (y<0.5f) y=0.5f; if (y>N+0.5f) y=N+0.5f; j0=(int)y; j1=j0+1;
		s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
		d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+
					 s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
	END_FOR
	set_bnd ( N, b, d );
}

void project ( int N, float * u, float * v, float * p, float * div )
{
	int i, j;

	FOR_EACH_CELL
		div[IX(i,j)] = -0.5f*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)])/N;
		p[IX(i,j)] = 0;
	END_FOR	
	set_bnd ( N, 0, div ); set_bnd ( N, 0, p );

	lin_solve ( N, 0, p, div, 1, 4 );

	FOR_EACH_CELL
		u[IX(i,j)] -= 0.5f*N*(p[IX(i+1,j)]-p[IX(i-1,j)]);
		v[IX(i,j)] -= 0.5f*N*(p[IX(i,j+1)]-p[IX(i,j-1)]);
	END_FOR
	set_bnd ( N, 1, u ); set_bnd ( N, 2, v );
}

void dens_step(int N, float *x, float *x0, float *u, float *v, float diff, float dt)
{
    add_source(N, x, x0, dt);
    SWAP(x0, x);
    diffuse(N, 0, x, x0, diff, dt);
    SWAP(x0, x);
    advect(N, 0, x, x0, u, v, dt);
}

void vel_step(int N, float *u, float *v, float *u0, float *v0, float visc, float dt)
{
    add_source(N, u, u0, dt);
    add_source(N, v, v0, dt);
    SWAP(u0, u);
    diffuse(N, 1, u, u0, visc, dt);
    SWAP(v0, v);
    diffuse(N, 2, v, v0, visc, dt);
    project(N, u, v, u0, v0);
    SWAP(u0, u);
    SWAP(v0, v);
    advect(N, 1, u, u0, u0, v0, dt);
    advect(N, 2, v, v0, u0, v0, dt);
    project(N, u, v, u0, v0);
}

/* global variables */

static int N;
static float dt, diff, visc;
static float force, source;
static int dvel;

static float *u, *v, *u_prev, *v_prev;
static float *dens, *dens_prev;

static void free_data(void)
{
    if (u)
        free(u);
    if (v)
        free(v);
    if (u_prev)
        free(u_prev);
    if (v_prev)
        free(v_prev);
    if (dens)
        free(dens);
    if (dens_prev)
        free(dens_prev);
}

static void clear_data(void)
{
    int i, size = (N + 2) * (N + 2);

    for (i = 0; i < size; i++)
    {
        u[i] = v[i] = u_prev[i] = v_prev[i] = dens[i] = dens_prev[i] = 0.0f;
    }
}

static int allocate_data(void)
{
    int size = (N + 2) * (N + 2);

    u = (float *)malloc(size * sizeof(float));
    v = (float *)malloc(size * sizeof(float));
    u_prev = (float *)malloc(size * sizeof(float));
    v_prev = (float *)malloc(size * sizeof(float));
    dens = (float *)malloc(size * sizeof(float));
    dens_prev = (float *)malloc(size * sizeof(float));

    if (!u || !v || !u_prev || !v_prev || !dens || !dens_prev)
    {
        fprintf(stderr, "cannot allocate data\n");
        return (0);
    }

    return (1);
}

static void get_from_UI(float *d, float *u, float *v)
{
    int i, j, size = (N + 2) * (N + 2);

    for (i = 0; i < size; i++)
    {
        u[i] = v[i] = d[i] = 0.0f;
    }

    i = N / 2;
    j = N / 2;

    if (i < 1 || i > N || j < 1 || j > N)
        return;

    u[IX(i, j)] = force * (1);
    v[IX(i, j)] = force * (1);
    d[IX(i, j)] = source;

    return;
}

static void printFrameMatrices(float *dens, float *u, float *v)
{
    for (int i = 0; i < 10; i++)
    {
        get_from_UI(dens_prev, u_prev, v_prev);
        vel_step(N, u, v, u_prev, v_prev, visc, dt);
        dens_step(N, dens, dens_prev, u, v, diff, dt);
        printf("Frame %d:\n", i);
        printf("Density Matrix:\t\t\t\t\tVelocity U Matrix:\t\t\t\t\tVelocity V Matrix:\n");
        for (int j = 0; j < N; j++)
        {
            printf("%f, ", dens[j * N]);
            printf("%f, ", dens[1 + j * N]);
            printf("%f, ", dens[2 + j * N]);
            printf("%f", dens[3 + j * N]);
            printf("\t\t");
            printf("%f, ", u[j * N]);
            printf("%f, ", u[1 + j * N]);
            printf("%f, ", u[2 + j * N]);
            printf("%f", u[3 + j * N]);
            if (i < 2 && j == 1)
            {
                printf("\t");
            }
            else
            {
                printf("\t\t");
            }
            if (j != 1)
            {
                printf("%f, ", v[j * N]);
                printf("%f, ", v[1 + j * N]);
                printf("%f, ", v[2 + j * N]);
                printf("%f", v[3 + j * N]);
            }
            else
            {
                printf("\t%f, ", v[j * N]);
                printf("%f, ", v[1 + j * N]);
                printf("%f, ", v[2 + j * N]);
                printf("%f", v[3 + j * N]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

int main(int argc, char **argv)
{
    // glutInit(&argc, argv);

    if (argc != 1 && argc != 7)
    {
        fprintf(stderr, "usage : %s N dt diff visc force source\n", argv[0]);
        fprintf(stderr, "where:\n");
        fprintf(stderr, "\t N      : grid resolution\n");
        fprintf(stderr, "\t dt     : time step\n");
        fprintf(stderr, "\t diff   : diffusion rate of the density\n");
        fprintf(stderr, "\t visc   : viscosity of the fluid\n");
        fprintf(stderr, "\t force  : scales the mouse movement that generate a force\n");
        fprintf(stderr, "\t source : amount of density that will be deposited\n");
        exit(1);
    }

    if (argc == 1)
    {
        N = 4;
        dt = 0.1f;
        diff = 0.0f;
        visc = 0.0f;
        force = 0.5f;
        source = 100.0f;
        fprintf(stderr, "Using defaults : N=%d dt=%g diff=%g visc=%g force = %g source=%g\n",
                N, dt, diff, visc, force, source);
    }
    else
    {
        N = atoi(argv[1]);
        dt = atof(argv[2]);
        diff = atof(argv[3]);
        visc = atof(argv[4]);
        force = atof(argv[5]);
        source = atof(argv[6]);
    }

    if (!allocate_data())
        exit(1);
    clear_data();

    printFrameMatrices(dens, u, v);

    clear_data();
    exit(0);
}