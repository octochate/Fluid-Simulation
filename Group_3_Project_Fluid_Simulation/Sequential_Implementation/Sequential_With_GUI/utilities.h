#pragma once
#ifndef _FLUIDSIM_UTILITIES
#define _FLUIDSIM_UTILITIES
#include <iostream>

#define MIN(x, y) (x > y) ? y : x
#define MAX(x, y) (x < y) ? y : x

class f3 {
public:
    float x;
    float y;
    float z;

    f3() {
        x = 0;
        y = 0;
        z = 0;
    }

    f3(float x_param, float y_param) {     // Constructor
        x = x_param;
        y = y_param;
        z = 0;
    }

    f3(float x_param, float y_param, float z_param) {
        x = x_param;
        y = y_param;
        z = z_param;
    }

    void cap(float maximumMagnitude) {
        float currentMagnitude = mag();
        if (currentMagnitude > maximumMagnitude) {
            float ratio = maximumMagnitude / currentMagnitude;
            x = x * ratio;
            y = y * ratio;
            z = z * ratio;
        }
    }

    static f3 abs(f3 obj) {
        return f3(std::abs(obj.x), std::abs(obj.y), std::abs(obj.z));
    }

    float magSquared() {
        return x * x + y * y + z * z;
    }

    float mag() {
        return sqrtf(x * x + y * y + z * z);
    }

    float dist(f3 c2) {
        return sqrtf(distSquared(c2));
    }

    float distSquared(f3 c2) {
        return (x - c2.x) * (x - c2.x) + (y - c2.y) * (y - c2.y) + (z - c2.z) * (z - c2.z);
    }

    //equality
    f3& operator=(const f3& c1);

    //addition
    friend f3 operator+(const f3& c1, const f3& c2);
    friend f3 operator+(const f3& c1, const float c2); //float
    friend f3 operator+(const float c2, const f3& c1);
    friend f3 operator+(const f3& c1, const double c2); //double
    friend f3 operator+(const double c2, const f3& c1);
    friend f3 operator+(const f3& c1, const int c2); //integer
    friend f3 operator+(const int c2, const f3& c1);

    //subtraction
    f3 operator-()const;
    friend f3 operator-(const f3& c1, const f3& c2);
    friend f3 operator-(const f3& c1, const float c2); //float
    friend f3 operator-(const float c2, const f3& c1);
    friend f3 operator-(const f3& c1, const double c2); //double
    friend f3 operator-(const double c2, const f3& c1);
    friend f3 operator-(const f3& c1, const int c2); //integer
    friend f3 operator-(const int c2, const f3& c1);

    //multiplication
    friend f3 operator*(const f3& c1, const f3& c2);
    friend f3 operator*(const f3& c1, const float c2); //float
    friend f3 operator*(const float c2, const f3& c1);
    friend f3 operator*(const f3& c1, const double c2); //double
    friend f3 operator*(const double c2, const f3& c1);
    friend f3 operator*(const f3& c1, const int c2); //integer
    friend f3 operator*(const int c2, const f3& c1);

    //division
    friend f3 operator/(const f3& c1, const f3& c2);
    friend f3 operator/(const f3& c1, const float c2); //float
    friend f3 operator/(const float c2, const f3& c1);
    friend f3 operator/(const f3& c1, const double c2); //double
    friend f3 operator/(const double c2, const f3& c1);
    friend f3 operator/(const f3& c1, const int c2); //integer
    friend f3 operator/(const int c2, const f3& c1);
};

std::ostream& operator<<(std::ostream& strm, const f3& a);

typedef struct GUI_force {
    f3 force;
    f3 position;
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

extern void force_enqueue(f3 force, f3 position, float impulseRadius);
extern GUI_force* force_dequeue();
extern int queueEmpty();
extern int queueSize();

extern void ink_enqueue(f3 force, f3 position, float impulseRadius);
extern GUI_force* ink_dequeue();
extern int ink_queueEmpty();
extern int ink_queueSize();
#endif