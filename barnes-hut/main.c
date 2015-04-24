// Richard Kirchofer
// barnes hut n-body simulation in n log n

// print
#include <stdio.h>
// malloc
#include <stdlib.h>

// 4 for quadtree, 8 for octree
#define CHILDREN 4
#define DIMENSIONS 2

struct node
{
	// array of 4 node pointers, one for each quadrant
	struct node* next[CHILDREN];
	// float coord[DIMENSIONS];
	float x;
	float y;
	float x_center;
	float y_center;
	// it may have a negative mass which means that it will repell instead of attract
	float mass;
	// width of the node (points don't have widths)
	float dimen;
};

// pass in how many child nodes for extensibility
void nullify_node( struct node* head )
// this function just nulls the values, not an init
{
	// nullify_node should not malloc memory because I will already have point structs
	// I pass in a point struct to add_point and it will get pointed to, no malloc
	for( int i = 0; i < CHILDREN; ++i)
	{
		head->next[i] = NULL;
	}
	head->x = 0.0;
	head->y = 0.0;
	// head.z = 0.0; // for 3d space
	head->x_center = 0.0;
	head->y_center = 0.0;
	head->mass = 0.0;
	head->dimen = 0.0;
}

// if next == nullptr and x&y == null then it is an external node
// an external node is an empty quadrant, it does not have any points
// nor does it have any other nodes
// if next == nullptr but x&y == values then it is a point and should have a mass
// if next != null and x&y != null then it is an internal node
// an internal node points to four other nodes
// x&y describe the position of the center of mass of the next nodes
void add_point( struct node* head, float dimen, struct node* point )
{
	// we can't yet add in our mass and average against the x,y
	// what if we already did that to the parent and we just need to make this a point
	if( head == NULL ) {
		// the quadrant is empty, we don't have external nodes, instead it's just null
		head = point;
	} else {
		// head is not null, then it should have an x,y and mass
		// it could still be a point though, we cant average ourselves just yet
		if( head->dimen != 0 ) {
			// ok, the width is nonzero, it's not a point
			// lets just add ourself to the register, "point was here"
			head->x = ((head->x * head->mass) + (point->x * point->mass)) / (head->mass + point->mass);
			head->y = ((head->y * head->mass) + (point->y * point->mass)) / (head->mass + point->mass);
			head->mass += point->mass;
			// we include ourself and now we need to recurse on the correct next pointer
			// compare point's x,y to find where it should go
			int quadrant = 0;
			quadrant ^= (point->x >= head->x_center);
			quadrant ^= (point->y >= head->y_center)<<1;
			add_point( head->next[quadrant], head->dimen / 2, point );
			/*
			// if to the right
			if( point->x >= x_center ) {
				// if down, right
				if( point->y <= y_center ) {
					add_point( head->next[1], point );
				} else {
				// if up, right
					add_point( head->next[0], point );
				}
			} else {
				// if to the left
				// if down, left
				if( point->y <= y_center ) {
					add_point( head->next[2], point );
				} else {
					// if up, left
					add_point( head->next[3], point );
				}
			}
			*/
		} else {
			// if the width is zero then it is a point
			// create new internal node and move existing point into new node
			// then move current point into new node
			struct node* old_point = head;
			head = (struct node*) malloc( sizeof( struct node ) );
			nullify_node( head );
			/*
			// average position and sum masses for head
			head->x = ((old_point->x * old_point->mass) + (point->x * point->mass)) / (old_point->mass + point->mass);
			head->y = ((old_point->y * old_point->mass) + (point->y * point->mass)) / (old_point->mass + point->mass);
			head->mass = old_point->mass + point->mass;
			*/
			// new node inherits dimension
			head->dimen = dimen / 2;
			// recurse on both old_point and point
			add_point( head, head->dimen, old_point );
			add_point( head, head->dimen, point );
		}
	//head = ( struct node* ) malloc( sizeof( struct node ) );
	}
}

int main(int argc, char* argv[])
{
	// create a pointer to a node, always keep track of this
	struct node* head;
	// dimen of head
	float dimen = 1.0;
	// test code for adding point
	struct node* point = (struct node*) malloc( sizeof( struct node ) );
	nullify_node( point );
	return 0;
}

/*
quadrants
0	|	1
____|____
	|
2	|	3
*//*	*//*
(0)(internal node)[1, 2, 3, null]
	/			|				\
(1)(point)		|				 \
		(2)(external node)		  \
					(3)(internal node)[13, 14, 15, 16]
						/			|			|		\
					(13)(point)		|			|		 \
								(14)(point)		|		  \
											(15)(point)	   \
														(16)(point)
*/

