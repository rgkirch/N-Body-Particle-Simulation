// Richard Kirchofer
// barnes hut n-body simulation in n log n

// has RAND_MAX
#include <cstdlib>
// atan2, sin, cos, fmod, pi
#include <math.h>
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
void nullify_node( struct node* current )
// this function just nulls the values, not an init
{
	// nullify_node should not malloc memory because I will already have point structs
	// I pass in a point struct to add_point and it will get pointed to, no malloc
	for( int i = 0; i < CHILDREN; ++i)
	{
		current->next[i] = nullptr;
	}
	current->x = 0.0;
	current->y = 0.0;
	// current.z = 0.0; // for 3d space
	current->x_center = 0.0;
	current->y_center = 0.0;
	current->mass = 0.0;
	current->dimen = 1.0;
	current->quadrant = 0;
}

// if next == nullptr and x&y == null then it is an external node
// an external node is an empty quadrant, it does not have any points
// nor does it have any other nodes
// if next == nullptr but x&y == values then it is a point and should have a mass
// if next != null and x&y != null then it is an internal node
// an internal node points to four other nodes
// x&y describe the position of the center of mass of the next nodes
void add_point( struct node* &previous, struct node* &current, struct node* &point )
{
	cout << "add point called" << endl;
	// we can't yet add in our mass and average against the x,y
	// what if we already did that to the parent and we just need to make this a point
	if( current == nullptr ) {
		cout << "addpt current == null" << endl;
		// the quadrant is empty, we don't have external nodes, instead it's just null
		current = point;
		cout << "addpt current address " << current << endl;
	} else {
		// current is not null, then it should have an x,y and mass
		// it could still be a point though, we cant average ourselves just yet
		if( current->dimen != 0 ) {
			cout << "addpt dimen != 0" << endl;
			// ok, the width is nonzero, it's not a point
			// lets just add ourself to the register, "point was here"
			current->x = ((current->x * current->mass) + (point->x * point->mass)) / (current->mass + point->mass);
			current->y = ((current->y * current->mass) + (point->y * point->mass)) / (current->mass + point->mass);
			current->mass += point->mass;
			// we include ourself and now we need to recurse on the correct next pointer
			// compare point's x,y to find where it should go
			point->quadrant = 0;
			point->quadrant ^= (point->x >= current->x_center);
			point->quadrant ^= (point->y >= current->y_center)<<1;
			add_point( current, current->next[point->quadrant], point );
		} else {
			cout << "addpt current dimen == 0" << endl;
			// if the width is zero then it is a point
			// create new internal node and move existing point into new node
			// then move current point into new node
			// reserve new memory for the internal node
			struct node* new_node = (struct node*) malloc( sizeof( struct node ) );
			cout << "addptr create new node" << endl;
			// zero out the new node
			nullify_node( new_node );
			// get quadrant from the old point
			new_node->quadrant = current->quadrant;
			current->quadrant = 0;
			// new node inherits dimension
			// nullptr from c++11
			new_node->dimen = previous->dimen / 2.0;
			previous->next[new_node->quadrant] = new_node;
			// if previous is null, then dimen should be 1.0 which is default
			// calculate the new center
			if( new_node->quadrant & 1 ) {
				new_node->x_center += new_node->dimen;
			} else {
				new_node->x_center -= new_node->dimen;
			}
			if( new_node->quadrant>>1 & 1 ) {
				new_node->y_center += new_node->dimen;
			} else {
				new_node->y_center -= new_node->dimen;
			}
			// recurse on both old_point and point
			// i'll need to recalculate the quadrant of the old point and the new point

			for( int iter = 0; iter < 4; ++iter ) {
				if( previous->next[iter] == new_node )
					cout << "previous->next[" << iter << "] == new node of " << new_node << endl;
			}
			for( int iter = 0; iter < 4; ++iter ) {
				if( new_node->next[iter] != nullptr )
					cout << "new node has non null" << endl;
			}
			add_point( previous, new_node, current);
			// TODO - if the x and y are the same, then there is infinite recursion
			add_point( previous, new_node, point );
			//add_point( previous, previous, new_node );
		}
	//current = ( struct node* ) malloc( sizeof( struct node ) );
	}
	cout << "addpt end" << endl;
}

// allocate memory for a point with some random x,y and return a pointer to it
struct node* random_point()
{
	struct node* current = (struct node*) malloc( sizeof( struct node ) );
	// default dimen is 1, point must have 0
	current->dimen = 0.0;
	current->mass= 1.0;
	// x,y between 0,1
	current->x = rand() / (float)RAND_MAX;
	current->y = rand() / (float)RAND_MAX;
	return current;
}

void draw( sf::RenderWindow &window, struct node* &current, sf::CircleShape &dot, sf::RectangleShape &box )
{
	cout << "draw called" << endl;
	if( current != nullptr ) {
		cout << "draw, current not nullptr" << endl;
		if( current->dimen > 0 ) {
			cout << "draw, current dimen > 0" << endl;
			// it's an internal node, let's draw a box
			box.setSize( sf::Vector2f( current->dimen * window.getSize().x, current->dimen * window.getSize().y ) );
			box.setPosition( current->x_center * window.getSize().x, current->y_center * window.getSize().y );
			cout << "draw box" << endl;
			window.draw( box );
			for( int iter = 0; iter < CHILDREN; ++iter ) {
				cout << "going to call on " << current->next[iter] << endl;
				draw( window, current->next[iter], dot, box );
			}
		} else {
			// it's a point, lets draw a dot
			cout << "draw dot" << endl;
			dot.setPosition( current->x * window.getSize().x, current->y * window.getSize().y);
			window.draw( dot );
		}
	}
}

void assign( struct node* &current, struct node* &point )
{
	current = point;
}

int main(int argc, char* argv[])
{
	// create a pointer to a node, always keep track of this
	struct node* head = (struct node*) malloc( sizeof( struct node ) );
	struct node* null = nullptr;
	nullify_node( head );
	struct node* dip = random_point();
	add_point( null, head, dip );
	dip = random_point();
	add_point( null, head, dip );
	dip = random_point();
	add_point( null, head, dip );

	/*
	random_pt = random_point();
	cout << "random point x y " << random_pt->x << " " << random_pt->y << endl;
	add_point( null, head, random_pt );
	cout << "random point quadrant " << random_pt->quadrant << endl;
	*/

	//head->next[random_pt->quadrant] = random_pt;

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
	box.setOutlineThickness( -8.0 );

	draw( window, head, dot, box );
	window.display();

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
		/*
		window.clear( sf::Color::White );
		// draw define elsewhere
		draw( window, head, dot, box );
		window.display();
		*/
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

