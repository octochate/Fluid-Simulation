#include <stdlib.h>
#include <stdio.h>
#include "utilities.h"


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