#include <assert.h>
// has RAND_MAX
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <vector>

#include <SFML/Graphics.hpp>

using namespace std;

#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 600
#define NSTEPS 500
#define SAVEFREQ 1

#define density	0.0005
#define mass	1.0
#define dt		0.0005
#define size 1.0

// double size = 0.0;

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

typedef struct 
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
} particle_t;

double read_timer( )
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

int find_option( int argc, char **argv, const char *option )
{
	for( int i = 1; i < argc; i++ )
		if( strcmp( argv[i], option ) == 0 )
			return i;
	return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
	int iplace = find_option( argc, argv, option );
	if( iplace >= 0 && iplace < argc-1 )
		return atoi( argv[iplace+1] );
	return default_value;
}

/*
void set_size( int n ) 
{
    size = sqrt( density * n );
}
*/

void init_particles( int n, particle_t *p )
{
	for( int i = 0; i < n; i++ ) 
	{
		p[i].x = 0.0;
		p[i].y = 0.0;
		p[i].vx = 0.0;
		p[i].vy = 0.0;
		p[i].ax = 0.0;
		p[i].ay = 0.0;
		// p[i].vx = rand()*1;
		// p[i].vy = rand()*1;
		p[i].x = rand() / (double)RAND_MAX;
		p[i].y = rand() / (double)RAND_MAX;
	}
}

void apply_force( particle_t &particle, particle_t &neighbor)
{
	double xdiff = neighbor.x - particle.x;
	double ydiff = neighbor.y - particle.y;
	double r2 = xdiff * xdiff + ydiff * ydiff;
	double force = 0;
	if (r2 != 0)
	{
		force = 1/r2;
	}
	double direction = atan2(xdiff, ydiff);
	particle.vx += sin(direction)*dt;
	particle.vy += cos(direction)*dt;
}

void move( particle_t &p )
{
	p.vx += p.ax * dt;
	p.vy += p.ay * dt;
	p.x  += p.vx * dt;
	p.y  += p.vy * dt;
}

void create_window(sf::RenderWindow &window)
{
	// create window
	window.create( sf::VideoMode( WINDOW_WIDTH, WINDOW_HEIGHT ), "SFML Visualizer!" );
	// window.setVerticalSyncEnabled( true );
	// window.setFramerateLimit( 15 );
}

int main( int argc, char **argv )
{	 
	int n = read_int( argc, argv, "-n", 1000 );
	//	set up MPI
	int n_proc, rank;
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	//	allocate memory for particles
	particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
	//	define a particle
	MPI_Datatype PARTICLE;
	MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
	MPI_Type_commit( &PARTICLE );
	//	set up the data partitioning across processors
	//	define the offsets
	int particle_per_proc = (n + n_proc - 1) / n_proc;
	int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
	for( int i = 0; i < n_proc+1; i++ )
		partition_offsets[i] = min( i * particle_per_proc, n );
	//	define the sizes
	int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
	for( int i = 0; i < n_proc; i++ )
		partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
	//	allocate storage for local partition
	int nlocal = partition_sizes[rank];
	particle_t *local = (particle_t*) malloc( nlocal * sizeof(particle_t) );
	//	initialize
	// size is a global
	// set_size( n );
	//	every process creates a pointer to a renderwindow
	sf::RenderWindow window;
	//	define a circle
	sf::CircleShape circle(4, 8);
	circle.setFillColor( sf::Color::Black );
	if( rank == 0 )
	{
		srand( time( NULL ) );
		//	init the particles
		init_particles( n, particles );
		//	only one process creates a window
		create_window(window);
	}
	//	distribute the particles
	MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
	//	simulate a number of time steps
	double simulation_time = read_timer( );
	// for( int step = 0; step < NSTEPS; step++ )
	while( 1 )
	{
		//	collect all global data locally
		MPI_Allgatherv( local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );
		if( rank == 0 )
		{
			window.clear( sf::Color::White );
			sf::CircleShape temp = circle;
			for( int i = 0; i < n; ++i )
			{
				temp.setPosition( particles[i].x / size * WINDOW_WIDTH, particles[i].y / size * WINDOW_HEIGHT );
				// cout << particles[i].x << " " << particles[i].y << endl;
				// temp.setPosition( .5 * WINDOW_WIDTH, .5 * WINDOW_HEIGHT );
				window.draw( temp );
			}
			window.display();
		}
		//	compute all forces
		for( int i = 0; i < nlocal; i++ )
		{
			local[i].ax = local[i].ay = 0;
			for (int j = 0; j < n; j++ )
				apply_force( local[i], particles[j]);
		}
	 
		//	move particles
		for( int i = 0; i < nlocal; i++ )
			move( local[i] );
	}
	simulation_time = read_timer( ) - simulation_time;
  
	if (rank == 0) {  
		printf( "n = %d, simulation time = %g seconds", n, simulation_time);
	}
  
	//	release resources
	free( partition_offsets );
	free( partition_sizes );
	free( local );
	free( particles );
	//	end mpi portion
	MPI_Finalize( );
	return 0;
}
