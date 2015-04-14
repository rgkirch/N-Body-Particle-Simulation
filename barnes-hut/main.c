// Richard Kirchofer
// barnes hut n-body simulation in n log n

// print
#include <stdio.h>
// malloc
#include <stdlib.h>

// 4 for quadtree, 8 for octree
#define CHILDREN 4

struct node
{
	// array of 4 node pointers, one for each quadrant
	struct node* next[CHILDREN];
	float x;
	float y;
	float mass;
	// width of the node (points don't have widths)
	float dimen;
};

// pass in how many child nodes for extensibility
void init_head( struct node* head, int children )
{
	for( int i = 0; i < children; ++i)
	{
		head->next[i] = NULL;
	}
	head->x = -1.0;
	head->y = -1.0;
	// head.z = -1; // for 3d space
	head->mass = 0.0;
}

// if next == nullptr and x&y == null then it is an external node
// an external node is an empty quadrant, it does not have any points
// nor does it have any other nodes
// if next == nullptr but x&y == values then it is a point and should have a mass
// if next != null and x&y != null then it is an internal node
// an internal node points to four other nodes
// x&y describe the position of the center of mass of the next nodes
void add_point( struct node* head, float x, float y, float mass )
{
	// not pointing to anything, either a point on an external node
	if( head->next == NULL )
	{
		// if there isn't mass, then let's say it's an external node (not a point) 
		// it may have a negative mass which means that it will repell instead of attract
		if( head->mass == 0 ) 
		{
			if( head->x != -1.0 ||
		}
	}
	//head = ( struct node* ) malloc( sizeof( struct node ) );
}

int main()
{
	struct node head;
	add_point( &head, 1.0, 1.0, 1.0 );
	return 0;
}
