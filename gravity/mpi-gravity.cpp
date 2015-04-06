#include <iostream>
#include <vector>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <SFML/Graphics.hpp>
#include "common-gravity.h"

using namespace std;

#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 600


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
	double size = set_size( n );
	// cout << "size:" << size << endl;
	//	every process creates a pointer to a renderwindow
	sf::RenderWindow window;
	//	define a circle
	sf::CircleShape circle(4, 8);
	circle.setFillColor( sf::Color::Black );
	if( rank == 0 )
	{
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
