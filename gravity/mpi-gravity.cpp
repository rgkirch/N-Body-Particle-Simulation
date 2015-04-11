// Richard Kirchofer

// for pi as a constant on older systems
#define _USE_MATH_DEFINES
// has RAND_MAX
#include <cstdlib>
// c++ input output
#include <iostream>
// atan2, sin, cos, fmod, pi
#include <math.h>
// everything
#include <mpi.h>
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

#include <SFML/Graphics.hpp>

using namespace std;

#define WINDOW_WIDTH 1920
#define WINDOW_HEIGHT 1080
//#define WINDOW_WIDTH 800
//#define WINDOW_HEIGHT 600


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
	double mass;
	sf::Color color;
} particle_t;

double read_timer();
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
void init_particles( int n, particle_t* p, int initialVelocity );
void apply_force( particle_t &particle, particle_t &neighbor, float cutoff, float timescale );
void move( particle_t &p, float timescale, float velFric );
int one_color( int color )
{
	int value = 0;
	int mod = color % 256;
	int div = color / 256;
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
	return value;
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
	// initial velocity random between 0 and this
	float initialVelocity = 10.0;
	// velocity drops to this fraction of itself every step
	float velFric = 0.99;
	// scale back how far they move, increases the quality of the calculations
	float timescale = 0.0005;
	// if they are farther away from eachother than this, no interaction besides attraction
	float cutoff = 0.00001;
	// number of particles
	int n = 100;
	// number of processors
	int n_proc;
	// processor id
	int rank;
	// how many steps to simulate
	int nSteps = 500;
	// should it save the frames or not
	int save_frames = 0;
	int rgb_array[3];
	rgb_array[0] = 0;
	rgb_array[1] = 0;
	rgb_array[2] = 0;
	/*
	vector< string > options;
	options.push_back("-n");
	options.push_back("--cutoff");
	options.push_back("--timescale");
	options.push_back("--velocityscale");
	options.push_back("--masstoradius");
	options.push_back("--defaultradius");
	options.push_back("--initialvelocity");
	options.push_back("--nsteps");
	vector< float > variables;
	variables.push_back(  );
	variables.push_back(  );
	variables.push_back(  );
	variables.push_back(  );
	variables.push_back(  );
	variables.push_back(  );
	variables.push_back(  );
	variables.push_back(  );
	*/
	/*
	for( int o = 1; o < argc; ++o )
	{
		if( strcmp( argv[o], "-n" ) == 0 )
			n = atoi(argv[o+1]);
		else
			n = 500;
		if( strcmp( argv[o], "--cutoff" ) == 0 )
			cutoff = argv[o+1];
		else
			cutoff = 0.00001;
		if( strcmp( argv[o], "--timescale" ) == 0 )
			timescale = argv[o+1];
		else
			timescale = 0.0005;
		if( strcmp( argv[o], "--velocityscale" ) == 0 )
			velocityscale = argv[o+1];
		else
			velocityscale = 0.99;
		if( strcmp( argv[o], "--masstoradius" ) == 0 )
			masstoradius = argv[o+1];
		else
			masstoradius = 1.0;
		if( strcmp( argv[o], "--defaultradius" ) == 0 )
			defaultradius = argv[o+1];
		else
			defaultradius = 4.0;
		if( strcmp( argv[o], "--initialvelocity" ) == 0 )
			initialvelocity = argv[o+1];
		else
			initialvelocity = 10;
		if( strcmp( argv[o], "--nsteps" ) == 0 )
			nsteps = argv[o+1];
		else
			nSteps = 50;
	}
	*/
	/*
	cout << "n " << n << endl;
	cout << "cutoff " << cutoff << endl;
	cout << "timescale " << timescale << endl;
	cout << "velFric " << velFric << endl;
	cout << "massToDrawRadius " << massToDrawRadius << endl;
	cout << "defaultRadius " << defaultRadius << endl;
	cout << "initialVelocity " << initialVelocity << endl;
	cout << "nSteps " << nSteps << endl;
	cout << endl;
	*/
	/*
	if( find_option( "-h" ) )
	{
		cout << "-h prints this help" << endl;
		cout << "-n number of particles in simulation" << endl;
	}
	*/

	//	set up MPI
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	//color_picker( rgb_array, rank, n_proc );
	//printf( "rank: %d, (%d, %d, %d)", rank, rgb_array[0], rgb_array[1], rgb_array[2] );
	//	allocate memory for particles
	// particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
	particle_t *particles = new particle_t[n];
	//	define a particle
	MPI_Datatype PARTICLE;
	MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
	MPI_Type_commit( &PARTICLE );
	//	set up the data partitioning across processors
	//	define the offsets
	int particle_per_proc = n / n_proc;
	// int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
	int* partition_offsets = new int[n_proc];
	for( int i = 0; i < n_proc; ++i )
		partition_offsets[i] = particle_per_proc * i;
	//	define the sizes
	// int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
	int* partition_sizes = new int[n_proc];
	for( int i = 0; i < n_proc; ++i )
		partition_sizes[i] = particle_per_proc;
	partition_sizes[n_proc-1] = n / n_proc + n % n_proc;
	//	allocate storage for local partition
	int nlocal = partition_sizes[rank];
	particle_t* local = new particle_t[nlocal];
	//	initialize
	//	every process creates a pointer to a renderwindow
	sf::RenderWindow window;
	// create a circle object, must be visible to all proc
	sf::CircleShape circle(1, 8);
	circle.setFillColor( sf::Color::Black );
	if( rank == 0 )
	{
		srand( time( NULL ) );
		//	init the particles
		init_particles( n, particles, initialVelocity );
		for( int i = 0; i < n; ++i )
		{
			particles[i].color = sf::Color::Blue;
		}
		//	only one process creates a window
		// function call defined above
		create_window(window);
	}
	//	distribute the particles
	MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
	//	simulate a number of time steps
	double simulation_time = read_timer( );
	for( int step = 0; step < nSteps; step++ )
	// while( 1 )
	{
		if( step == 0)
		{
			for( int i = 0; i < nlocal; ++i )
			{
				local[i].color = sf::Color::Red;
			}
		}
		//	collect all global data locally
		MPI_Allgatherv( local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );
		// cout << "P " << partition_offsets[n_proc-1] << endl;
		if( rank == 0 )
		{
			window.clear( sf::Color::White );
			for( int i = 0; i < n; ++i )
			{
				// V = 4/3*pi*r**3
				// pow((3/4)V/M_PI, 1/3)
				// pow(3.0/4.0/M_PI*particles[i].mass,1.0/3.0)
				circle.setRadius(particles[i].mass);
				// circle.setRadius(pow(3.0/4.0/M_PI*particles[i].mass,1.0/3.0));
				circle.setFillColor(particles[i].color);
				// fmod is floating point modulus
				// don't do modulus TODO
				circle.setPosition( fmod( particles[i].x * WINDOW_WIDTH, WINDOW_WIDTH), fmod(particles[i].y * WINDOW_HEIGHT, WINDOW_HEIGHT ) );
				// cout << particles[i].x << " " << particles[i].y << endl;
				// circle.setPosition( .5 * WINDOW_WIDTH, .5 * WINDOW_HEIGHT );
				window.draw( circle );
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
		}
		//	compute all forces
		for( int i = 0; i < nlocal; i++ )
		{
			local[i].ax = local[i].ay = 0;
			for (int j = 0; j < n; j++ )
				apply_force( local[i], particles[j], cutoff, timescale );
		}
	 
		//	move particles
		for( int i = 0; i < nlocal; i++ )
		{
			move( local[i], timescale, velFric );
		}
	}
	simulation_time = read_timer( ) - simulation_time;
  
	if (rank == 0) {  
		printf( "n = %d, simulation time = %g seconds", n, simulation_time);
	}
  
	//	release resources
	// free( partition_offsets );
	delete partition_offsets;
	// free( partition_sizes );
	delete partition_sizes;
	// free( local );
	delete local;
	// free( particles );
	delete particles;
	//	end mpi portion
	MPI_Finalize( );
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
		p[i].vx = rand() / (double)RAND_MAX * initialVelocity;
		p[i].vy = rand() / (double)RAND_MAX * initialVelocity;
		p[i].x = rand() / (double)RAND_MAX;
		p[i].y = rand() / (double)RAND_MAX;
		p[i].color = sf::Color::Black;
	}
}

void apply_force( particle_t &particle, particle_t &neighbor, float cutoff, float timescale)
{
	double dx = neighbor.x - particle.x;
	double dy = neighbor.y - particle.y;
	double r2 = dx * dx + dy * dy;
	if( r2 == 0 ) {
		return;
	} else if( r2 < cutoff ) {
		//double xdiff = (double)abs(particle.vx - neighbor.vx);
		//double ydiff = (double)abs(particle.vy - neighbor.vy);
		double xsum = (particle.vx + neighbor.vx);
		double ysum = (particle.vy + neighbor.vy);
		double xavg = xsum / 2;
		double yavg = ysum / 2;
		particle.vx = xavg;
		particle.vy = yavg;
	} else {
		double force = 0.0;
		force = 1/r2;
		double direction = atan2(dx, dy);
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
}

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
