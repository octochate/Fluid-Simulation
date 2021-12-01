#include <stdlib.h>
#include <stdio.h>
#include "utilities.h"

std::ostream& operator<<(std::ostream& strm, const f3& a) {
	return strm << "f2(" << a.x << ", " << a.y << ", " << a.z << ")";
}

f3& f3::operator=(const f3& c1) { x = c1.x; y = c1.y; z = c1.z;  return *this; } //f2

////////////////
// addition
// ////////////
f3 operator+(const f3& c1, const f3& c2) { return f3(c1.x + c2.x, c1.y + c2.y, c1.z + c2.z); } //f2
f3 operator+(const f3& c1, const float c2) { return f3(c1.x + c2, c1.y + c2, c1.z + c2); }  //float
f3 operator+(const float c2, const f3& c1) { return f3(c1.x + c2, c1.y + c2, c1.z + c2); }
f3 operator+(const f3& c1, const double c2) { return f3(c1.x + (float)c2, c1.y + (float)c2, c1.z + (float)c2); } //double
f3 operator+(const double c2, const f3& c1) { return f3(c1.x + (float)c2, c1.y + (float)c2, c1.z + (float)c2); }
f3 operator+(const f3& c1, const int c2) { return f3(c1.x + c2, c1.y + c2, c1.z + c2); } //int
f3 operator+(const int c2, const f3& c1) { return f3(c1.x + c2, c1.y + c2, c1.z + c2); }

////////////////
// subtraction
// ////////////
f3 operator-(const f3& c1, const f3& c2) { return f3(c1.x - c2.x, c1.y - c2.y, c1.z - c2.z); } //f2
f3 f3::operator-()const { return f3(-x, -y, -z); }
f3 operator-(const f3& c1, const float c2) { return f3(c1.x - c2, c1.y - c2, c1.z - c2); } //float
f3 operator-(const float c2, const f3& c1) { return f3(c1.x - c2, c1.y - c2, c1.z - c2); }
f3 operator-(const f3& c1, const double c2) { return f3(c1.x - (float)c2, c1.y - (float)c2, c1.z - (float)c2); } //double
f3 operator-(const double c2, const f3& c1) { return f3(c1.x - (float)c2, c1.y - (float)c2, c1.z - (float)c2); }
f3 operator-(const f3& c1, const int c2) { return f3(c1.x - c2, c1.y - c2, c1.z - c2); } //int
f3 operator-(const int c2, const f3& c1) { return f3(c1.x - c2, c1.y - c2, c1.z - c2); }

////////////////
// multiplication
// ////////////
f3 operator*(const f3& c1, const f3& c2) { return f3(c1.x * c2.x, c1.y * c2.y, c1.z * c2.z); } //f2
f3 operator*(const f3& c1, const float c2) { return f3(c1.x * c2, c1.y * c2, c1.z * c2); } //float
f3 operator*(const float c2, const f3& c1) { return f3(c1.x * c2, c1.y * c2, c1.z * c2); }
f3 operator*(const f3& c1, const double c2) { return f3(c1.x * (float)c2, c1.y * (float)c2, c1.z * (float)c2); } //double
f3 operator*(const double c2, const f3& c1) { return f3(c1.x * (float)c2, c1.y * (float)c2, c1.z * (float)c2); }
f3 operator*(const f3& c1, const int c2) { return f3(c1.x * c2, c1.y * c2, c1.z * c2); } //int
f3 operator*(const int c2, const f3& c1) { return f3(c1.x * c2, c1.y * c2, c1.z * c2); }

////////////////
// division
// ////////////
f3 operator/(const f3& c1, const f3& c2) { return f3(c1.x / c2.x, c1.y / c2.y, c1.z / c2.z); } //f2
f3 operator/(const f3& c1, const float c2) { return f3(c1.x / c2, c1.y / c2, c1.z / c2); } //float
f3 operator/(const float c2, const f3& c1) { return f3(c1.x / c2, c1.y / c2, c1.z / c2); }
f3 operator/(const f3& c1, const double c2) { return f3(c1.x / (float)c2, c1.y / (float)c2, c1.z / (float)c2); } //double
f3 operator/(const double c2, const f3& c1) { return f3(c1.x / (float)c2, c1.y / (float)c2, c1.z / (float)c2); }
f3 operator/(const f3& c1, const int c2) { return f3(c1.x / c2, c1.y / c2, c1.z / c2); } //int
f3 operator/(const int c2, const f3& c1) { return f3(c1.x / c2, c1.y / c2, c1.z / c2); }


queue forceQueue = { NULL, NULL, 0 };
queue inkQueue = { NULL, NULL, 0 };


extern void force_enqueue(f3 force, f3 position, float impulseRadius) {
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

extern void ink_enqueue(f3 force, f3 position, float impulseRadius) {
	node* newNode = (node*)malloc(sizeof(node));
	newNode->force = (GUI_force*)malloc(sizeof(GUI_force));
	newNode->force->force = force;
	newNode->force->position = position;
	newNode->force->impulseRadius = impulseRadius;
	newNode->next = NULL;
	newNode->previous = NULL;

	if (ink_queueSize() == 0) {
		inkQueue.first = newNode;
		inkQueue.last = newNode;
	}
	else {
		inkQueue.last->next = newNode;
		newNode->previous = inkQueue.last;
		inkQueue.last = newNode;
	}
	inkQueue.size++;
}

extern GUI_force* ink_dequeue() {
	if (ink_queueSize() == 0) {
		return NULL;
	}
	GUI_force* force;
	if (ink_queueSize() == 1) {
		//get last element from queue
		force = inkQueue.last->force;

		//remove last element from queue, which is also the first element of queue
		free(inkQueue.last);
		inkQueue.first = NULL;
		inkQueue.last = NULL;
	}
	else {
		//get last element from queue
		force = inkQueue.last->force;

		//remove last element from queue
		node* toDelete = inkQueue.last;
		inkQueue.last = inkQueue.last->previous;
		inkQueue.last->next = NULL;
		free(toDelete);
	}
	inkQueue.size--;
	return force;
}

extern int ink_queueEmpty() {
	return ink_queueSize() == 0;
}

extern int ink_queueSize() {
	return inkQueue.size;
}