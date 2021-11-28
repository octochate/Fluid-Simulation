#include <stdlib.h>
#include <stdio.h>
#include "utilities.h"

std::ostream& operator<<(std::ostream& strm, const f2& a) {
	return strm << "f2(" << a.x << ", " << a.y << ")";
}

f2& f2::operator=(const f2& c1) { x = c1.x; y = c1.y; return *this; } //f2

////////////////
// addition
// ////////////
f2 operator+(const f2& c1, const f2& c2) { return f2(c1.x + c2.x, c1.y + c2.y); } //f2
f2 operator+(const f2& c1, const float c2) { return f2(c1.x + c2, c1.y + c2); }  //float
f2 operator+(const float c2, const f2& c1) { return f2(c1.x + c2, c1.y + c2); }
f2 operator+(const f2& c1, const double c2) { return f2(c1.x + (float)c2, c1.y + (float)c2); } //double
f2 operator+(const double c2, const f2& c1) { return f2(c1.x + (float)c2, c1.y + (float)c2); }
f2 operator+(const f2& c1, const int c2) { return f2(c1.x + c2, c1.y + c2); } //int
f2 operator+(const int c2, const f2& c1) { return f2(c1.x + c2, c1.y + c2); }

////////////////
// subtraction
// ////////////
f2 operator-(const f2& c1, const f2& c2) { return f2(c1.x - c2.x, c1.y - c2.y); } //f2
f2 f2::operator-()const { return f2(-x, -y); }
f2 operator-(const f2& c1, const float c2) { return f2(c1.x - c2, c1.y - c2); } //float
f2 operator-(const float c2, const f2& c1) { return f2(c1.x - c2, c1.y - c2); }
f2 operator-(const f2& c1, const double c2) { return f2(c1.x - (float)c2, c1.y - (float)c2); } //double
f2 operator-(const double c2, const f2& c1) { return f2(c1.x - (float)c2, c1.y - (float)c2); }
f2 operator-(const f2& c1, const int c2) { return f2(c1.x - c2, c1.y - c2); } //int
f2 operator-(const int c2, const f2& c1) { return f2(c1.x - c2, c1.y - c2); }

////////////////
// multiplication
// ////////////
f2 operator*(const f2& c1, const f2& c2) { return f2(c1.x * c2.x, c1.y * c2.y); } //f2
f2 operator*(const f2& c1, const float c2) { return f2(c1.x * c2, c1.y * c2); } //float
f2 operator*(const float c2, const f2& c1) { return f2(c1.x * c2, c1.y * c2); }
f2 operator*(const f2& c1, const double c2) { return f2(c1.x * (float)c2, c1.y * (float)c2); } //double
f2 operator*(const double c2, const f2& c1) { return f2(c1.x * (float)c2, c1.y * (float)c2); }
f2 operator*(const f2& c1, const int c2) { return f2(c1.x * c2, c1.y * c2); } //int
f2 operator*(const int c2, const f2& c1) { return f2(c1.x * c2, c1.y * c2); }

////////////////
// division
// ////////////
f2 operator/(const f2& c1, const f2& c2) { return f2(c1.x / c2.x, c1.y / c2.y); } //f2
f2 operator/(const f2& c1, const float c2) { return f2(c1.x / c2, c1.y / c2); } //float
f2 operator/(const float c2, const f2& c1) { return f2(c1.x / c2, c1.y / c2); }
f2 operator/(const f2& c1, const double c2) { return f2(c1.x / (float)c2, c1.y / (float)c2); } //double
f2 operator/(const double c2, const f2& c1) { return f2(c1.x / (float)c2, c1.y / (float)c2); }
f2 operator/(const f2& c1, const int c2) { return f2(c1.x / c2, c1.y / c2); } //int
f2 operator/(const int c2, const f2& c1) { return f2(c1.x / c2, c1.y / c2); }


queue forceQueue = { NULL, NULL, 0 };

extern void force_enqueue(f2 force, f2 position, float impulseRadius) {
	node* newNode = (node*)malloc(sizeof(node));
	newNode->force = (GUI_force*)malloc(sizeof(GUI_force));
	newNode->force->force = force;
	newNode->force->position = position;
	newNode->force->impulseRadius = impulseRadius;
	newNode->next = NULL;
	newNode->previous = NULL;

	if (queueSize() == 0) {
		forceQueue.first = newNode;
		forceQueue.last = newNode;
	}
	else {
		forceQueue.last->next = newNode;
		newNode->previous = forceQueue.last;
		forceQueue.last = newNode;
	}
	forceQueue.size++;
}

extern GUI_force* force_dequeue() {
	if (queueSize() == 0) {
		return NULL;
	}
	GUI_force* force;
	if (queueSize() == 1) {
		//get last element from queue
		force = forceQueue.last->force;

		//remove last element from queue, which is also the first element of queue
		free(forceQueue.last);
		forceQueue.first = NULL;
		forceQueue.last = NULL;
	}
	else {
		//get last element from queue
		force = forceQueue.last->force;

		//remove last element from queue
		node* toDelete = forceQueue.last;
		forceQueue.last = forceQueue.last->previous;
		forceQueue.last->next = NULL;
		free(toDelete);
	}
	forceQueue.size--;
	return force;
}

extern int queueEmpty() {
	return queueSize() == 0;
}

extern int queueSize() {
	return forceQueue.size;
}