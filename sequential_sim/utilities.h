#pragma once
#ifndef _FLUIDSIM_UTILITIES
#define _FLUIDSIM_UTILITIES
#include <iostream>

#define MIN(x, y) (x > y) ? y : x
#define MAX(x, y) (x < y) ? y : x

class f2 {
public:
    float x;
    float y;

    f2() {
        x = 0;
        y = 0;
    }

    f2(float x_param, float y_param) {     // Constructor
        x = x_param;
        y = y_param;
    }

    void cap(float min, float max) {
        x = MAX(min, MIN(x, max));
        y = MAX(min, MIN(y, max));
    }

    //equality
    f2& operator=(const f2& c1);

    //addition
    friend f2 operator+(const f2& c1, const f2& c2);
    friend f2 operator+(const f2& c1, const float c2); //float
    friend f2 operator+(const float c2, const f2& c1);
    friend f2 operator+(const f2& c1, const double c2); //double
    friend f2 operator+(const double c2, const f2& c1);
    friend f2 operator+(const f2& c1, const int c2); //integer
    friend f2 operator+(const int c2, const f2& c1);

    //subtraction
    f2 operator-()const;
    friend f2 operator-(const f2& c1, const f2& c2);
    friend f2 operator-(const f2& c1, const float c2); //float
    friend f2 operator-(const float c2, const f2& c1);
    friend f2 operator-(const f2& c1, const double c2); //double
    friend f2 operator-(const double c2, const f2& c1);
    friend f2 operator-(const f2& c1, const int c2); //integer
    friend f2 operator-(const int c2, const f2& c1);

    //multiplication
    friend f2 operator*(const f2& c1, const f2& c2);
    friend f2 operator*(const f2& c1, const float c2); //float
    friend f2 operator*(const float c2, const f2& c1);
    friend f2 operator*(const f2& c1, const double c2); //double
    friend f2 operator*(const double c2, const f2& c1);
    friend f2 operator*(const f2& c1, const int c2); //integer
    friend f2 operator*(const int c2, const f2& c1);

    //division
    friend f2 operator/(const f2& c1, const f2& c2);
    friend f2 operator/(const f2& c1, const float c2); //float
    friend f2 operator/(const float c2, const f2& c1);
    friend f2 operator/(const f2& c1, const double c2); //double
    friend f2 operator/(const double c2, const f2& c1);
    friend f2 operator/(const f2& c1, const int c2); //integer
    friend f2 operator/(const int c2, const f2& c1);
};

std::ostream& operator<<(std::ostream& strm, const f2& a);

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