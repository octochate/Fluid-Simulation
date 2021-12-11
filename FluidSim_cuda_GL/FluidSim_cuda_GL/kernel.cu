#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda_runtime_api.h>
#include <time.h>
#include "freeglut.h"
#include <mutex>

/* macros */
#define IX(i,j) ((i)+(N+2)*(j))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define MIN(x, y) (x > y) ? y : x
#define MAX(x, y) (x < y) ? y : x

static int allocate_data(void);

__global__ void add_source(float* x, float* s, float dt, int size)
{
	int idx = (blockIdx.x * blockDim.x + threadIdx.x);

	if (idx < size) {
		x[idx] += dt * s[idx];
	}
}

__device__ void set_bnd(int N, int b, float* x, int index, int elementsPerThread)
{
	int i = index + 1;
	if (i > N + 1) {
		return;
	}

	int size = (N + 2) * (N + 2);

	while (i < (index + elementsPerThread) && i <= N + 1) {
		x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
		x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
		x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
		x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
		i++;
	}

	__syncthreads();

	if (index == 0) {
		x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
		x[IX(0, N + 1)] = 0.5f * (x[IX(1, N + 1)] + x[IX(0, N)]);
		x[IX(N + 1, 0)] = 0.5f * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
		x[IX(N + 1, N + 1)] = 0.5f * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
	}
}

__global__ void lin_solve(int N, int b, float* x, float* x0, float a, float c, int elementPerThread)
{
	int i, j, k, idxNew;
	int idx = (blockIdx.x * blockDim.x + threadIdx.x) * elementPerThread;

	if (idx >= N * N) {
		return;
	}

	for (k = 0; k < 20; k++)
	{
		idxNew = idx;
		while (idxNew < idx + elementPerThread && idxNew < N * N) {
			j = idxNew / N + 1;
			i = idxNew % N + 1;
			x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] + x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
			idxNew++;
		}
		__syncthreads();
		set_bnd(N, b, x, idx, elementPerThread);
	}
}

__global__ void advect(int N, int b, float* d, float* d0, float* u, float* v, float dt, int elementPerThread)
{
	int idx = (blockIdx.x * blockDim.x + threadIdx.x) * elementPerThread;

	if (idx >= N * N) {
		return;
	}

	int idxNew = idx;
	while (idxNew < idx + elementPerThread && idxNew < N * N) {
		int j = idxNew / N + 1;
		int i = idxNew % N + 1;

		int i0, j0, i1, j1;
		float x, y, s0, t0, s1, t1, dt0;

		dt0 = dt * N;
		x = i - dt0 * u[IX(i, j)]; y = j - dt0 * v[IX(i, j)];
		if (x < 0.5f) x = 0.5f; if (x > N + 0.5f) x = N + 0.5f; i0 = (int)x; i1 = i0 + 1;
		if (y < 0.5f) y = 0.5f; if (y > N + 0.5f) y = N + 0.5f; j0 = (int)y; j1 = j0 + 1;
		s1 = x - i0; s0 = 1 - s1; t1 = y - j0; t0 = 1 - t1;

		d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
			s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);

		idxNew++;
	}

	__syncthreads();
	set_bnd(N, b, d, idx, elementPerThread);
}

__global__ void project1(int N, float* u, float* v, float* p, float* div, int elementPerThread)
{
	int idx = (blockIdx.x * blockDim.x + threadIdx.x) * elementPerThread;

	if (idx >= N * N) {
		return;
	}

	int idxNew = idx;
	while (idxNew < idx + elementPerThread && idxNew < N * N) {
		int j = idxNew / N + 1;
		int i = idxNew % N + 1;

		div[IX(i, j)] = -0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
		p[IX(i, j)] = 0;
		idxNew++;
	}

	__syncthreads();
	set_bnd(N, 0, div, idx, elementPerThread);
	set_bnd(N, 0, p, idx, elementPerThread);
}

__global__ void project3(int N, float* u, float* v, float* p, int elementPerThread)
{
	int idx = (blockIdx.x * blockDim.x + threadIdx.x) * elementPerThread;

	if (idx >= N * N) {
		return;
	}

	int idxNew = idx;
	while (idxNew < idx + elementPerThread && idxNew < N * N) {
		int j = idxNew / N + 1;
		int i = idxNew % N + 1;

		u[IX(i, j)] -= 0.5f * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
		v[IX(i, j)] -= 0.5f * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);

		idxNew++;
	}
	__syncthreads();
	set_bnd(N, 1, u, idx, elementPerThread);
	set_bnd(N, 2, v, idx, elementPerThread);
}


/* global variables */
cudaStream_t stream1 = NULL, stream2 = NULL, stream3 = NULL;

volatile static int N, size;
volatile static float dt, diff, visc, force, source;
volatile static int display_velocity = 0, cuda_streams = 0, shoot_liquid = 0;

float numIntervals = 100;
float dtMax = 5;
float dtMin = 0.01;
float dt_del = (dtMax - dtMin) / numIntervals;
float viscMax = 1;
float viscMin = 0;
float visc_del = (viscMax - viscMin) / numIntervals;
float diffMax = 1.;
float diffMin = 0;
float diff_del = (diffMax - diffMin) / numIntervals;
float fluidAmountMax = 1000.;
float fluidAmountMin = 1.;
float fluidAmount_del = (fluidAmountMax - fluidAmountMin) / numIntervals;
float forceAmountMax = 5;
float forceAmountMin = 0.1;
float forceAmount_del = (forceAmountMax - forceAmountMin) / numIntervals;

// Parallelize computation :: TODO update values to match
int num_threads_source;
int num_blocks_source;
int num_threads;
int elementsPerThread;

// CPU variables
static float* u, * v, * dens;
static float* u_userInput, * v_userInput, * dens_userInput;

volatile int flag = 0;
std::mutex * guimutexPtr;

// GPU variables
float* u_cuda, * v_cuda, * u_prev_cuda, * v_prev_cuda;
float* p_cuda, * div_cuda;
float* dens_cuda, * dens_prev_cuda;

static int win_id;
static int win_x, win_y;
static int mouse_down[3];
static int omx, omy, mx, my;
static float xtext, ytext;


/*
  ----------------------------------------------------------------------
   free/clear/allocate simulation data
  ----------------------------------------------------------------------
*/


static void free_data(void)
{
	if (u_userInput) free(u_userInput);
	if (v_userInput) free(v_userInput);
	if (dens_userInput) free(dens_userInput);

	if (u) free(u);
	if (v) free(v);
	if (dens) free(dens);

	if (u_cuda) cudaFree(u_cuda);
	if (v_cuda) cudaFree(v_cuda);
	if (u_prev_cuda) cudaFree(u_prev_cuda);
	if (v_prev_cuda) cudaFree(v_prev_cuda);
	if (p_cuda) cudaFree(p_cuda);
	if (div_cuda) cudaFree(div_cuda);
	if (dens_cuda) cudaFree(dens_cuda);
	if (dens_prev_cuda) cudaFree(dens_prev_cuda);
}

static void clear_data(void)
{
	cudaMemset(u_cuda, 0, sizeof(float) * size);
	cudaMemset(v_cuda, 0, sizeof(float) * size);
	cudaMemset(p_cuda, 0, sizeof(float) * size);
	cudaMemset(div_cuda, 0, sizeof(float) * size);
	cudaMemset(dens_cuda, 0, sizeof(float) * size);
	cudaMemset(u_prev_cuda, 0, sizeof(float) * size);
	cudaMemset(v_prev_cuda, 0, sizeof(float) * size);
	cudaMemset(dens_prev_cuda, 0, sizeof(float) * size);

	for (int i = 0; i < size; i++) {
		u_userInput[i] = v_userInput[i] = dens_userInput[i] = 0.0f;
	}
}

static void destroy_streams(void) {
	// Destroy cuda streams
	if (stream1 == stream2 && stream2 == stream3) {
		cudaStreamDestroy(stream1);
	}
	else {
		cudaStreamDestroy(stream1);
		cudaStreamDestroy(stream2);
		cudaStreamDestroy(stream3);
	}

	stream1 = NULL;
	stream2 = NULL;
	stream3 = NULL;
}

static void create_streams(void) {
	if (stream1 != NULL) {
		return;
	}

	if (cuda_streams) {
		cudaStreamCreate(&stream1);
		cudaStreamCreate(&stream2);
		cudaStreamCreate(&stream3);
	}
	else {
		cudaStreamCreate(&stream1);
		stream2 = stream1;
		stream3 = stream1;
	}
}

static void set_numThreads(int newValue) {
	if (newValue < 1) {
		newValue = 1;
	}
	if (newValue > 1024) {
		newValue = 1024;
	}

	num_threads = newValue;
	elementsPerThread = (N * N + 1) / num_threads;

}

static int set_gridSize(int newValue) {
	if (newValue < 32) {
		newValue = 32;
	}
	if (newValue > 1024) {
		newValue = 1024;
	}

	if (N == newValue) {
		return 1;
	}

	//cudaDeviceReset();
	cudaDeviceSynchronize();
	const std::lock_guard<std::mutex> lock(*guimutexPtr);
	free_data();
	N = newValue;
	size = (N + 2) * (N + 2);
	num_threads_source = (N + 2);
	num_blocks_source = (N + 2);
	set_numThreads(num_threads);
	return allocate_data();
}

static int allocate_data(void)
{
	// Allocate space for device copies
	u = (float*)malloc(size * sizeof(float));
	v = (float*)malloc(size * sizeof(float));
	dens = (float*)malloc(size * sizeof(float));

	u_userInput = (float*)malloc(size * sizeof(float));
	v_userInput = (float*)malloc(size * sizeof(float));
	dens_userInput = (float*)malloc(size * sizeof(float));

	// gpu copies
	cudaMalloc(&u_cuda, sizeof(float) * size);
	cudaMalloc(&v_cuda, sizeof(float) * size);
	cudaMalloc(&p_cuda, sizeof(float) * size);
	cudaMalloc(&div_cuda, sizeof(float) * size);
	cudaMalloc(&u_prev_cuda, sizeof(float) * size);
	cudaMalloc(&v_prev_cuda, sizeof(float) * size);
	cudaMalloc(&dens_cuda, sizeof(float) * size);
	cudaMalloc(&dens_prev_cuda, sizeof(float) * size);


	if (!u_userInput || !v_userInput || !dens_userInput || !dens || !u || !v || !u_cuda || !v_cuda || !u_prev_cuda || !v_prev_cuda || !dens_cuda || !dens_prev_cuda) {
		fprintf(stderr, "cannot allocate data\n");
		return (0);
	}

	return (1);
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

	h = 1.0f / N;

	cudaMemcpy(u, u_cuda, sizeof(float) * size, cudaMemcpyDeviceToHost);
	cudaMemcpy(v, v_cuda, sizeof(float) * size, cudaMemcpyDeviceToHost);

	glColor3f(1.0f, 1.0f, 1.0f);
	glLineWidth(1.0f);

	glBegin(GL_LINES);

	for (i = 1; i <= N; i++) {
		x = (i - 0.5f) * h;
		for (j = 1; j <= N; j++) {
			y = (j - 0.5f) * h;

			glVertex2f(x, y);
			glVertex2f(x + u[IX(i, j)], y + v[IX(i, j)]);
		}
	}

	glEnd();
}

static void draw_density(void)
{
	int i, j;
	float x, y, h, d00, d01, d10, d11;

	h = 1.0f / N;
	cudaMemcpy(dens, dens_cuda, sizeof(float) * size, cudaMemcpyDeviceToHost);

	glBegin(GL_QUADS);

	for (i = 0; i <= N; i++) {
		x = (i - 0.5f) * h;
		for (j = 0; j <= N; j++) {
			y = (j - 0.5f) * h;

			d00 = dens[IX(i, j)];
			d01 = dens[IX(i, j + 1)];
			d10 = dens[IX(i + 1, j)];
			d11 = dens[IX(i + 1, j + 1)];

			glColor3f(d00, d00, d00 * 0); glVertex2f(x, y);
			glColor3f(d10, d10, d10 * 0); glVertex2f(x + h, y);
			glColor3f(d11, d11, d11 * 0); glVertex2f(x + h, y + h);
			glColor3f(d01, d01, d01 * 0); glVertex2f(x, y + h);
		}
	}

	glEnd();
}

/*
  ----------------------------------------------------------------------
   relates mouse movements to forces sources
  ----------------------------------------------------------------------
*/

static void get_from_UI(float* d, float* u, float* v)
{
	cudaMemcpy(u_prev_cuda, u_userInput, sizeof(float) * size, cudaMemcpyHostToDevice);
	cudaMemcpy(v_prev_cuda, v_userInput, sizeof(float) * size, cudaMemcpyHostToDevice);;
	cudaMemcpy(dens_prev_cuda, dens_userInput, sizeof(float) * size, cudaMemcpyHostToDevice);
	return;
}

static void get_from_UI_CPU(void)
{
	int i, j;

	if (shoot_liquid) {
		i = (N + 2) / 5;
		j = (N + 2) / 5;

		u_userInput[IX(i, j)] = force;
		v_userInput[IX(i, j)] = force;
		dens_userInput[IX(i, j)] = source;
	}

	if (!mouse_down[0] && !mouse_down[2]) return;

	i = (int)((mx / (float)win_x) * N + 1);
	j = (int)(((win_y - my) / (float)win_y) * N + 1);

	if (i<1 || i>N || j<1 || j>N) return;

	if (mouse_down[0]) {
		u_userInput[IX(i, j)] += force * (mx - omx);
		v_userInput[IX(i, j)] += force * (omy - my);
	}

	if (mouse_down[2]) {
		dens_userInput[IX(i - 1, j)] = source;
		dens_userInput[IX(i + 1, j)] = source;
		dens_userInput[IX(i, j + 1)] = source;
		dens_userInput[IX(i, j - 1)] = source;
		dens_userInput[IX(i, j)] = source;
	}

	omx = mx;
	omy = my;

	return;
}

/*
  ----------------------------------------------------------------------
   GLUT callback routines
  ----------------------------------------------------------------------
*/
static void Sliders()
{
	const unsigned char a[50] = "Controls";
	const unsigned char b[50] = "Time Step    [a]";
	const unsigned char c[50] = "Viscosity    [s]";
	const unsigned char d[50] = "Diffussion   [d]";
	const unsigned char e[50] = "Fluid Amount [f]";
	const unsigned char f[50] = "Force Amount [g]";
	const unsigned char g[50] = "# threads: ";
	const unsigned char h[50] = "grid size: ";

	unsigned char numThreads[50];
	strcpy((char*) numThreads, std::to_string(num_threads).c_str());
	const unsigned char* numThreads_const = (const unsigned char*)numThreads;

	unsigned char gridSize[50];
	strcpy((char*)gridSize, std::to_string(N).c_str());
	const unsigned char* gridSize_const = (const unsigned char*)gridSize;

	const unsigned char cudaStreams[50] = "CUDA streams enabled";
	const unsigned char dispenser[50] = "Jet enabled";

	const unsigned char* aPtr = a;
	const unsigned char* bPtr = b;
	const unsigned char* cPtr = c;
	const unsigned char* dPtr = d;
	const unsigned char* ePtr = e;
	const unsigned char* fPtr = f;
	const unsigned char* gPtr = g;
	const unsigned char* hPtr = h;
	const unsigned char* streamsPtr = cudaStreams;
	const unsigned char* jetPtr = dispenser;

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
	glRasterPos2f(0.78125 + 0.05, 0.7109375);
	glutBitmapString(GLUT_BITMAP_8_BY_13, gPtr);
	glRasterPos2f(0.78125 + 0.05 + 0.1, 0.7109375);
	glutBitmapString(GLUT_BITMAP_8_BY_13, numThreads_const);
	glRasterPos2f(0.78125 + 0.05, 0.6718750 + 0.025);
	glutBitmapString(GLUT_BITMAP_8_BY_13, hPtr);
	glRasterPos2f(0.78125 + 0.05 + 0.1, 0.6718750 + 0.025);
	glutBitmapString(GLUT_BITMAP_8_BY_13, gridSize_const);
	
	if (cuda_streams) {
		glRasterPos2f(0.78125, 0.6718750);
		glutBitmapString(GLUT_BITMAP_8_BY_13, streamsPtr);
	}
	if (shoot_liquid) {
		glRasterPos2f(0.78125, 0.6718750 - 0.025);
		glutBitmapString(GLUT_BITMAP_8_BY_13, jetPtr);
	}

	glRasterPos2f(0., 0.);

	glBegin(GL_LINES);
	glColor3f(1.0, 1.0, 1.0);

	// Draw slider boxes.
	for (int i = 0; i < 5; i++)
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

	// Compute dynamic slider fill.
	float dtSliderEnd = ((dt / dtMax) * 0.15625) + sliderStart;
	float viscSliderEnd = ((visc / viscMax) * 0.15625) + sliderStart;
	float diffSliderEnd = ((diff / diffMax) * 0.15625) + sliderStart;
	float fluidAmountSliderEnd = ((source / fluidAmountMax) * 0.15625) + sliderStart;
	float forceAmountSliderEnd = ((force / forceAmountMax) * 0.15625) + sliderStart;

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
		if (i <= diffSliderEnd)
		{
			heightTop = 1. - (38. + 2. * 20.) / 512.;
			heightBottom = 1. - (49. + 2. * 20.) / 512.;
			glVertex2d(i, heightTop);
			glVertex2d(i, heightBottom);
		}
		if (i <= fluidAmountSliderEnd)
		{
			heightTop = 1. - (38. + 3. * 20.) / 512.;
			heightBottom = 1. - (49. + 3. * 20.) / 512.;
			glVertex2d(i, heightTop);
			glVertex2d(i, heightBottom);
		}
		if (i <= forceAmountSliderEnd)
		{
			heightTop = 1. - (38. + 4. * 20.) / 512.;
			heightBottom = 1. - (49. + 4. * 20.) / 512.;
			glVertex2d(i, heightTop);
			glVertex2d(i, heightBottom);
		}
	}
	glEnd();
}

static void key_func(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'c':
	case 'C':
		clear_data();
		break;

	case 'v':
	case 'V':
		display_velocity = !display_velocity;
		break;

	case 'a':
		dt = MIN(dtMax, dt+ dt_del);
		printf("dt is now %f\n", dt);
		break;

	case 'A':
		dt = MAX(dtMin, dt - dt_del);
		printf("dt is now %f\n", dt);
		break;
	case 's':
		visc = MIN(viscMax, visc + visc_del);
		printf("visc is now %f\n", visc);
		break;

	case 'S':
		visc = MAX(viscMin, visc - visc_del);
		printf("visc is now %f\n", visc);
		break;
	case 'd':
		diff = MIN(diffMax, diff + diff_del);
		printf("diff is now %f\n", diff);
		break;

	case 'D':
		diff = MAX(diffMin, diff - diff_del);
		printf("diff is now %f\n", diff);
		break;
	case 'f':
		source = MIN(fluidAmountMax, source + fluidAmount_del);
		printf("fluidAmount is now %f\n", source);
		break;

	case 'F':
		source = MAX(fluidAmountMin, source - fluidAmount_del);
		printf("fluidAmount is now %f\n", source);
		break;
	case 'g':
		force = MIN(forceAmountMax, force + forceAmount_del);
		printf("forceAmount is now %f\n", force);
		break;

	case 'G':
		force = MAX(forceAmountMin, force - forceAmount_del);
		printf("forceAmount is now %f\n", force);
		break;

	case 'e':
	case 'E':
		destroy_streams();
		cuda_streams = !cuda_streams;
		create_streams();
		break;
	case 'w':
	case 'W':
		shoot_liquid = !shoot_liquid;
		break;
	case '1':
		set_numThreads(num_threads * 2);
		break;
	case '2':
		set_numThreads(num_threads / 2);
		break;
	case '3':
		if (!set_gridSize(N + 100)) exit(1);
		clear_data();
		break;
	case '4':
		if (!set_gridSize(N - 100)) exit(1);
		clear_data();
		break;
	}

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

static void reshape_func(int width, int height)
{
	glutSetWindow(win_id);
	glutReshapeWindow(width, height);

	win_x = width;
	win_y = height;
}

static void idle_func(void)
{
	// have copy of host data
	get_from_UI(dens_prev_cuda, u_prev_cuda, v_prev_cuda);
	cudaDeviceSynchronize();

	////// Velocity timestep parallelization
	add_source << < num_blocks_source, num_threads_source >> > (u_cuda, u_prev_cuda, dt, size);
	add_source << < num_blocks_source, num_threads_source >> > (v_cuda, v_prev_cuda, dt, size);
	add_source << < num_blocks_source, num_threads_source >> > (dens_cuda, dens_prev_cuda, dt, size);

	for (int i = 0; i < size; i++) {
		u_userInput[i] = v_userInput[i] = dens_userInput[i] = 0.0f;
	}

	SWAP(u_prev_cuda, u_cuda);
	SWAP(v_prev_cuda, v_cuda);
	SWAP(dens_prev_cuda, dens_cuda);
	cudaDeviceSynchronize();

	

	//// diffuse
	float a = dt * visc * N * N;
	lin_solve << < 1, num_threads, 0, stream1 >> > (N, 1, u_cuda, u_prev_cuda, a, 1 + 4 * a, elementsPerThread); // diffuse u_cuda
	lin_solve << < 1, num_threads, 0, stream2 >> > (N, 2, v_cuda, v_prev_cuda, a, 1 + 4 * a, elementsPerThread); // diffuse v_cuda

	a = dt * diff * N * N;
	lin_solve << < 1, num_threads, 0, stream3 >> > (N, 0, dens_cuda, dens_prev_cuda, a, 1 + 4 * a, elementsPerThread); // diffuse dens_cuda
	
	cudaDeviceSynchronize();

	// projection step (no swapping beforehand)
	project1 << < 1, num_threads >> > (N, u_cuda, v_cuda, p_cuda, div_cuda, elementsPerThread);
	lin_solve << < 1, num_threads >> > (N, 0, p_cuda, div_cuda, 1, 4, elementsPerThread);
	project3 << < 1, num_threads >> > (N, u_cuda, v_cuda, p_cuda, elementsPerThread);
	SWAP(u_prev_cuda, u_cuda);
	SWAP(v_prev_cuda, v_cuda);
	cudaDeviceSynchronize();

	advect << < 1, num_threads, 0, stream1 >> > (N, 1, u_cuda, u_prev_cuda, u_prev_cuda, v_prev_cuda, dt, elementsPerThread);
	advect << < 1, num_threads, 0, stream2 >> > (N, 2, v_cuda, v_prev_cuda, u_prev_cuda, v_prev_cuda, dt, elementsPerThread);
	cudaDeviceSynchronize();

	// projection step (no swapping beforehand)
	project1 << < 1, num_threads >> > (N, u_cuda, v_cuda, p_cuda, div_cuda, elementsPerThread);
	lin_solve << < 1, num_threads >> > (N, 0, p_cuda, div_cuda, 1, 4, elementsPerThread);
	project3 << < 1, num_threads >> > (N, u_cuda, v_cuda, p_cuda, elementsPerThread);
	SWAP(dens_prev_cuda, dens_cuda);
	cudaDeviceSynchronize();

	// Density timestep parallelization
	advect << < 1, num_threads >> > (N, 0, dens_cuda, dens_prev_cuda, u_cuda, v_cuda, dt, elementsPerThread);

	glutSetWindow(win_id);
	glutPostRedisplay();
}

static void display_func(void)
{
	pre_display();

	if (display_velocity) draw_velocity();
	else		draw_density();

	Sliders();
	post_display();
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
	// glutReshapeFunc ( reshape_func );
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
		//N = 128;
		dt = 1.0f;
		diff = 0.0f;
		visc = 0.0f;
		force = 1.0f;
		source = 100.0f;
		fprintf(stderr, "Using defaults : N=%d dt=%g diff=%g visc=%g force = %g source=%g\n",
			N, dt, diff, visc, force, source);
	}
	else {
		//N = atoi(argv[1]);
		dt = atof(argv[2]);
		diff = atof(argv[3]);
		visc = atof(argv[4]);
		force = atof(argv[5]);
		source = atof(argv[6]);
	}

	printf("\n\nHow to use this demo:\n\n");
	printf("\t Add densities with the right mouse button\n");
	printf("\t Add velocities with the left mouse button and dragging the mouse\n");
	printf("\t Toggle density/velocity display with the 'v' key\n");
	printf("\t Clear the simulation by pressing the 'c' key\n");
	printf("\t Quit by pressing the 'q' key\n");

	set_numThreads(1024);
	display_velocity = 0;

	std::mutex guimutex;
	guimutexPtr = &guimutex;

	if (!set_gridSize(300)) exit(1);
	clear_data();
	create_streams();

	
	// A mutex ensures orderly access to std::cout from multiple threads.
	std::thread t1([&guimutex]() {
		while (!flag) {
			{
				const std::lock_guard<std::mutex> lock(guimutex);
				get_from_UI_CPU();
			}
			std::this_thread::sleep_for(std::chrono::milliseconds(25));
		}
		});


	win_x = 800;
	win_y = 800;
	open_glut_window();

	glutMainLoop();

	//stop GUI thread
	flag = 1;
	t1.join();

	// Free Device space
	cudaFree(u_cuda);
	cudaFree(v_cuda);
	cudaFree(u_prev_cuda);
	cudaFree(dens_cuda);
	cudaFree(dens_prev_cuda);
	cudaFree(v_prev_cuda);

	destroy_streams();

	exit(0);
}