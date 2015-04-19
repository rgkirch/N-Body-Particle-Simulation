// Richard Kirchofer
// barnes hut n-body simulation in n log n

// print
#include <stdio.h>
#include <iostream>
// malloc
#include <stdlib.h>
#include <SFML/Graphics.hpp>

using namespace std;

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
	unsigned int quadrant;
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
	head->dimen = 1.0;
	head->quadrant = 0;
}

// if next == nullptr and x&y == null then it is an external node
// an external node is an empty quadrant, it does not have any points
// nor does it have any other nodes
// if next == nullptr but x&y == values then it is a point and should have a mass
// if next != null and x&y != null then it is an internal node
// an internal node points to four other nodes
// x&y describe the position of the center of mass of the next nodes
void add_point( struct node* previous, struct node* current, struct node* point )
{
	// we can't yet add in our mass and average against the x,y
	// what if we already did that to the parent and we just need to make this a point
	if( current == NULL ) {
		// the quadrant is empty, we don't have external nodes, instead it's just null
		current = point;
	} else {
		// current is not null, then it should have an x,y and mass
		// it could still be a point though, we cant average ourselves just yet
		if( current->dimen != 0 ) {
			// ok, the width is nonzero, it's not a point
			// lets just add ourself to the register, "point was here"
			current->x = ((current->x * current->mass) + (point->x * point->mass)) / (current->mass + point->mass);
			current->y = ((current->y * current->mass) + (point->y * point->mass)) / (current->mass + point->mass);
			current->mass += point->mass;
			// we include ourself and now we need to recurse on the correct next pointer
			// compare point's x,y to find where it should go
			point->quadrant ^= (point->x >= current->x_center);
			point->quadrant ^= (point->y >= current->y_center)<<1;
			add_point( current, current->next[point->quadrant], point );
		} else {
			// if the width is zero then it is a point
			// create new internal node and move existing point into new node
			// then move current point into new node
			// save the place of the current point
			struct node* old_point = current;
			// reserve new memory for the internal node
			current = (struct node*) malloc( sizeof( struct node ) );
			// zero out the new node
			nullify_node( current );
			// get quadrant from the old point
			current->quadrant = old_point->quadrant;
			// new node inherits dimension
			// nullptr from c++11
			if( previous != nullptr ) {
				current->dimen = previous->dimen / 2.0;
			}
			// if previous is null, then dimen should be 1.0 which is default
			// calculate the new center
			if( current->quadrant & 1 ) {
				current->x_center += current->dimen;
			} else {
				current->x_center -= current->dimen;
			}
			if( current->quadrant<<1 & 1 ) {
				current->y_center += current->dimen;
			} else {
				current->y_center -= current->dimen;
			}
			// recurse on both old_point and point
			// i'll need to recalculate the quadrant of the old point and the new point
			add_point( previous, current, old_point );
			// TODO - if the x and y are the same, then there is infinite recursion
			add_point( previous, current, point );
		}
	//current = ( struct node* ) malloc( sizeof( struct node ) );
	}
}

void draw( sf::RenderWindow &window, struct node* current, sf::CircleShape &dot, sf::RectangleShape &box )
{
	if( current != nullptr ) {
		if( current->dimen != 0 ) {
			// it's an internal node, let's draw a box
			box.setSize( sf::Vector2f( current->dimen, current->dimen ) );
			box.setPosition( current->x_center, current->y_center );
			window.draw( box );
			for( int iter = 0; iter < CHILDREN; ++iter ) {
				draw( window, current->next[iter], dot, box );
			}
		} else {
			// it's a point, lets draw a dot
			dot.setPosition( current->x, current->y);
			window.draw( dot );
		}
	}
}

int main(int argc, char* argv[])
{
	// create a pointer to a node, always keep track of this
	struct node* head = (struct node*) malloc( sizeof( struct node ) );
	nullify_node( head );

	// these values represent the default windowed dimensions (not fullscreen)
	int window_width = 800;
	int window_height = 600;

	sf::RenderWindow window;
	// limit the refresh rate to 15 fps, only works in fullscreen, uses sleeps
	//window.setFramerateLimit( 15 );
	window.create( sf::VideoMode( window_width, window_height ), "barnes_hut" );
	window.clear( sf::Color::White );

	sf::CircleShape dot;
	sf::RectangleShape box;

	dot.setRadius( 4 );
	dot.setFillColor( sf::Color::Black );

	box.setSize( sf::Vector2f( window.getSize().x, window.getSize().y ) );
	box.setPosition( sf::Vector2f( window.getSize().x / 2.0, window.getSize().y / 2.0 ) );
	cout << window.getSize().x << " " << window.getSize().y << endl;
	box.setFillColor( sf::Color::Transparent );
	box.setOutlineColor( sf::Color::Black );
	// a negative thickness means that it expands inwards towards the center of the shape
	box.setOutlineThickness( -20.0 );

	while( window.isOpen() )
	{
		sf::Event event;
		while( window.pollEvent( event ) )
		{
			if( event.type == sf::Event::Closed ) {
				// cout << "closing window" << endl;
				window.close();
			}
			if( event.type == sf::Event::LostFocus ) {
				// cout << "focus lost" << endl;
			}
			if( event.type == sf::Event::GainedFocus ) {
				// cout << "focus gained" << endl;
			}
			if( event.type == sf::Event::MouseMoved) {
				// cout << sf::Mouse::getPosition().x << " " << sf::Mouse::getPosition().y << endl;
			}
			if( event.type == sf::Event::MouseButtonPressed) {
				// cout << "mouse button pressed" << endl;
			}
			if( event.type == sf::Event::Resized ) {
				// cout << "resized to " << event.size.width << " by " << event.size.height << endl;
				window.setView( sf::View( sf::FloatRect( 0, 0, event.size.width, event.size.height ) ) );
				window.clear( sf::Color::White );
			}
		}
		window.clear( sf::Color::White );
		// draw define elsewhere
		draw( window, head, dot, box );
		window.display();
	}
	return 0;
}

/*
quadrant
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

