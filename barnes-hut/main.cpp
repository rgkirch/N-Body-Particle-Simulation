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

void draw( sf::RenderWindow &window, struct node* head, &dot, &box )
{
	dot.setRadius( 4 );
	dot.setFillColor( sf::Color::Black );

	box.setFillColor( sf::Color::Black );
	box.setPosition();
	box.setSize( sf::Vector2f( 8, 8 ) );

	if( bar.getPosition().x >= window.getSize().x )
	{
		window.clear( sf::Color::White );
		bar.setPosition(0, 0.5 * window.getSize().y);
	}
	window.draw( bar );
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

	box.setSize( sf::Vector2f( 8, 8 ) );
	box.setFillColor( sf::Color::Black );

	while( window.isOpen() )
	{
		sf::Event event;
		while( window.pollEvent( event ) )
		{
			if( event.type == sf::Event::Closed ) {
				cout << "closing window" << endl;
				window.close();
			}
			if( event.type == sf::Event::LostFocus ) {
				cout << "focus lost" << endl;
			}
			if( event.type == sf::Event::GainedFocus ) {
				cout << "focus gained" << endl;
			}
			if( event.type == sf::Event::MouseMoved) {
				cout << sf::Mouse::getPosition().x << " " << sf::Mouse::getPosition().y << endl;
			}
			if( event.type == sf::Event::MouseButtonPressed) {
				cout << "mouse button pressed" << endl;
			}
			if( event.type == sf::Event::Resized ) {
				cout << "resized to " << event.size.width << " by " << event.size.height << endl;
				window.setView( sf::View( sf::FloatRect( 0, 0, event.size.width, event.size.height ) ) );
				window.clear( sf::Color::White );
			}
		}
		draw( window, head, dot, box );
		window.display();
	}
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

