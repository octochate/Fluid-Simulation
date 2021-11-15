#pragma once
#ifndef _FLUIDSIM_UTILITIES
#define _FLUIDSIM_UTILITIES

typedef struct f2 {
    float x;
    float y;
}f2;

typedef struct GUI_force {
    f2 force;
    f2 position;
    float impulseRadius;
}GUI_force;

typedef struct node {
    struct node* next;
    struct node* previous;
    GUI_force* force;
}node;

typedef struct {
    node* first;
    node* last;
    int size;
} queue;

extern queue forceQueue;

extern void force_enqueue(f2 force, f2 position, float impulseRadius);
extern GUI_force* force_dequeue();
extern int queueEmpty();
extern int queueSize();
#endif