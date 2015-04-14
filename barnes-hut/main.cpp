// Richard Kirchofer
// barnes hut n-body simulation in n log n

#include <iostream>
#include <stdio.h>

using namespace std;

// if next == nullptr and x&y == null then it is an external node
// an external node is an empty quadrant, it does not have any points
// nor does it have any other nodes
// if next == nullptr but x&y == values then it is a point and should have a mass
// if next != null and x&y != null then it is an internal node
// an internal node points to four other nodes
// x&y describe the position of the center of mass of the next nodes
struct node
{
	// array of 4 node pointers, one for each quadrant
	struct node* next[4];
	int x;
	int y;
	int mass;
};

int main()
{
	struct node root { .next=nullptr, .x=-1, .y=-1, .mass=0 }
	return 0;
}
