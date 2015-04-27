// Richard Kirchofer

// for pi as a constant on older systems
#define _USE_MATH_DEFINES
// has RAND_MAX
#include <cstdlib>
// c++ input output
#include <iostream>
// atan2, sin, cos, fmod, pi
#include <math.h>
// ostringstream for converting int to string for frame_name
#include <sstream>
// i like printf
#include <stdio.h>
// strcmp
#include <string.h>
#include <sys/time.h>
#include <time.h>
// have a list of the available options
#include <vector>
#include "omp.h"

#include <SFML/Graphics.hpp>

using namespace std;

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

typedef struct 
{
	float x;
	float y;
	float vx;
	float vy;
	float ax;
	float ay;
	float mass;
	sf::Color color;
} particle_t;

float read_timer();
void init_particles( int n, particle_t* p, int initialVelocity );
void apply_force( particle_t &particle, particle_t &neighbor, float cutoff, float timescale );
void move( particle_t &p, float timescale, float velFric );
int one_color( int color )
{
	int value = 0;
	int mod = color % 256;
	int div = color / 256;
	printf( "%d \n", div );
	if( div == 2 ) {
		// if increasing
		value = mod;
	} else if( div == 5 ) {
		// if decreasing
		value = 255 - mod;
	} else if( div == 3 || div == 4 ) {
		value = 255;
	} // div == 0 or 1, value = 0; div != 6 because max color = 6*256
}

int* color_picker( int* value, int rank, int n_proc )
{
	if( rank >= n_proc )
	{
		cerr << "Invalid rank " << rank << ". Must be less than " << n_proc << endl;
	}
	int range = 6*256;
	int color = range / n_proc * rank;
	value[0] = one_color( color );
	value[1] = one_color( color + 2 * 256 );
	value[2] = one_color( color + 4 * 256 );
	printf( "rank: %d, (%d, %d, %d)", rank, value[0], value[1], value[2] );
	return value;
}

int main( int argc, char *argv[] )
{	 
	// these values represent the default windowed dimensions (not fullscreen)
	int window_width = 1920;
	int window_height = 1080;
	// initial velocity random between 0 and this
	float initialVelocity = 00.0;
	// velocity drops to this fraction of itself every step
	float velFric = 0.999;
	// scale back how far they move, increases the quality of the calculations
	float timescale = 0.0005;
	// if they are farther away from eachother than this, no interaction besides attraction
	float cutoff = 0.000001;
	// number of particles
	int n = 1000;
	// number of processors
	int n_proc = 0;
	//omp_set_num_threads(4);
	#pragma omp parallel
	if( omp_get_thread_num() == 0 )
		n_proc = omp_get_num_threads();
	printf( "n_proc: %d\n", n_proc );
	// how many steps to simulate
	int nSteps = 100;
	// should it save the frames or not
	int save_frames = 0;
	int rgb_array[3];
	rgb_array[0] = 0;
	rgb_array[1] = 0;
	rgb_array[2] = 0;

	//color_picker( rgb_array, rank, n_proc );
	//	allocate memory for particles
	particle_t *particles = new particle_t[n];
	// offsets
	//int particle_per_proc = ceil(n / n_proc);
	int particle_per_proc = (n + n_proc - 1) / n_proc;
	int* partition_offsets = new int[n_proc+1];
    for( int i = 0; i < n_proc+1; i++ )
        partition_offsets[i] = min( i * particle_per_proc, n );
	// sizes
	int* partition_sizes = new int[n_proc];
    for( int i = 0; i < n_proc; i++ )
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
	//	initialize
	sf::RenderWindow window;
	// create a circle object
	sf::CircleShape circle(2, 8);
	circle.setFillColor( sf::Color::Black );
	srand( time( NULL ) );
	//	init the particles
	init_particles( n, particles, initialVelocity );
	window.create( sf::VideoMode( window_width, window_height ), "gravity" );
	window.clear( sf::Color::White );
	float simulation_time = read_timer( );
	int step = 0;

	for( int i = 0; i < omp_get_num_threads(); ++i )
	{
		//printf( "%d\n", partition_sizes[i] );
	}

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
				//cout << "resized to " << event.size.width << " by " << event.size.height << endl;
				printf( "resized to %d by %d\n", event.size.width, event.size.height );
				window.setView( sf::View( sf::FloatRect( 0, 0, event.size.width, event.size.height ) ) );
				window.clear( sf::Color::White );
			}
		}
		window.clear( sf::Color::White );
		//#pragma omp parallel
		for( int i = 0; i < n; ++i )
		{
			// V = 4/3*pi*r**3
			// pow((3/4)V/M_PI, 1/3)
			// pow(3.0/4.0/M_PI*particles[i].mass,1.0/3.0)
			//circle.setRadius(particles[i].mass);
			// circle.setRadius(pow(3.0/4.0/M_PI*particles[i].mass,1.0/3.0));
			//circle.setFillColor(particles[i].color);
			// fmod is floating point modulus
			// don't do modulus TODO
			circle.setPosition( particles[i].x * window.getSize().x, particles[i].y * window.getSize().y );
			// cout << particles[i].x << " " << particles[i].y << endl;
			// circle.setPosition( .5 * WINDOW_WIDTH, .5 * WINDOW_HEIGHT );
			//float before = read_timer();
			//#pragma omp critical
			window.draw( circle );
			//printf( "%g\n", read_timer() - before );
		}
		if( save_frames )
		{
			sf::Image frame = window.capture();
			ostringstream stream;
			stream << "frame" << step << ".png";
			string frame_name = stream.str();
			frame.saveToFile(frame_name);
		}
		window.display();
		//	compute all forces
		#pragma omp parallel
		for( int i = partition_offsets[omp_get_thread_num()]; i < partition_offsets[omp_get_thread_num()+1]; ++i )
		{
			particles[i].ax = particles[i].ay = 0;
			for (int j = 0; j < n; j++ )
				apply_force( particles[i], particles[j], cutoff, timescale );
		}
	 
		//	move particles
		#pragma omp parallel
		for( int i = partition_offsets[omp_get_thread_num()]; i < partition_offsets[omp_get_thread_num()+1]; ++i )
		{
			move( particles[i], timescale, velFric );
		}
		++step;
		//printf( "%d\n", step );
		/*
		if( step >= nSteps )
			window.close();
		*/
	}
	printf( "simulation time: %g\n", read_timer() - simulation_time );
  
	//	release resources
	// free( particles );
	delete particles;
	return 0;
}

void init_particles( int n, particle_t* p, int initialVelocity )
{
	for( int i = 0; i < n; i++ ) 
	{
		p[i].x = 0.0;
		p[i].y = 0.0;
		p[i].vx = 0.0;
		p[i].vy = 0.0;
		p[i].ax = 0.0;
		p[i].ay = 0.0;
		p[i].mass = 3.0;
		// p[i].displaySize = p[i].mass * 4.0;
		p[i].vx = (rand() / (float)RAND_MAX * initialVelocity);
		p[i].vy = (rand() / (float)RAND_MAX * initialVelocity);
		p[i].x = rand() / (float)RAND_MAX;
		p[i].y = rand() / (float)RAND_MAX;
		p[i].color = sf::Color::Black;
	}
}

void apply_force( particle_t &particle, particle_t &neighbor, float cutoff, float timescale)
{
	float dx = neighbor.x - particle.x;
	float dy = neighbor.y - particle.y;
	float r2 = dx * dx + dy * dy;
	if( r2 == 0 ) {
		// printf( "%p %p\n", &particle, &neighbor );
	} else if( r2 < cutoff ) {
		//float xdiff = (float)abs(particle.vx - neighbor.vx);
		//float ydiff = (float)abs(particle.vy - neighbor.vy);
		float xsum = (particle.vx + neighbor.vx);
		float ysum = (particle.vy + neighbor.vy);
		float xavg = xsum / 2;
		float yavg = ysum / 2;
		particle.vx = xavg;
		particle.vy = yavg;
	} else {
		float force = 0.0;
		force = 1/r2;
		float direction = atan2(dx, dy);
		particle.vx += sin(direction)*timescale;
		particle.vy += cos(direction)*timescale;
	}
}

void move( particle_t &p, float timescale, float velFric )
{
	p.vx += p.ax * timescale;
	p.vy += p.ay * timescale;
	p.vx *= velFric;
	p.vy *= velFric;
	p.x  += p.vx * timescale;
	p.y  += p.vy * timescale;
	if( p.x < 0 )
		p.x += 1;
	if( p.x > 1 )
		p.x -= 1;
	if( p.y < 0 )
		p.y += 1;
	if( p.y > 1 )
		p.y -= 1;
}

float read_timer( )
{
	static bool initialized = false;
	static struct timeval start;
	struct timeval end;
	if( !initialized )
	{
		gettimeofday( &start, NULL );
		initialized = true;
	}
	gettimeofday( &end, NULL );
	return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}
