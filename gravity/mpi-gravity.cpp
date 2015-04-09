// Richard Kirchofer

// has RAND_MAX
#include <cstdlib>
// c++ input output
#include <iostream>
// atan2, sin, cos, fmod
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

// #define WINDOW_WIDTH 1920
// #define WINDOW_HEIGHT 1080
#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 600

float massToDrawRadius;
float defaultRadius;
float initialVelocity;
float velFric;
float timescale;
float cutoff;
int n;
int n_proc;
int rank;
int nSteps;

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
	double displaySize;
	sf::Color color;
} particle_t;

double read_timer();
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
void init_particles( int n, particle_t* p );
void apply_force( particle_t &particle, particle_t &neighbor);
void move( particle_t &p );

void create_window(sf::RenderWindow &window)
{
	// create window
	window.create( sf::VideoMode( WINDOW_WIDTH, WINDOW_HEIGHT ), "SFML Visualizer!" );
	// window.setVerticalSyncEnabled( true );
	// window.setFramerateLimit( 15 );
}

int main( int argc, char **argv )
{	 
	vector< string > options;
	options.push_back("-n");
	options.push_back("--cutoff");
	options.push_back("--timescale");
	options.push_back("--velocityscale");
	options.push_back("--masstoradius");
	options.push_back("--defaultradius");
	options.push_back("--initialvelocity");
	options.push_back("--nsteps");


	for( int v = 1; v < argc; ++v )
	{
		for( int o = 0; o < options.size(); ++o )
		{
			if( strcmp( argv[v], options[o]) == 0 )
			{
				
			}
		}
	}
	n = read_int( argc, argv, "-n", 1000 );
	cutoff = read_int( argc, argv, "--cutoff", 0.00001 );
	timescale = read_int( argc, argv, "--timescale", 0.0005 );
	velFric = read_int( argc, argv, "--velocityscale", 0.99 );
	massToDrawRadius = read_int(argc, argv, "--masstoradius", 1.0);
	defaultRadius = read_int(argc, argv, "--defaultradius", 4.0);
	initialVelocity = read_int(argc, argv, "--initialvelocity", 10);
	nSteps = read_int(argc, argv, "--nsteps", 50);
	cout << "n " << n << endl;
	cout << "cutoff " << cutoff << endl;
	cout << "timescale " << timescale << endl;
	cout << "velFric " << velFric << endl;
	cout << "massToDrawRadius " << massToDrawRadius << endl;
	cout << "defaultRadius " << defaultRadius << endl;
	cout << "initialVelocity " << initialVelocity << endl;
	cout << "nSteps " << nSteps << endl;
	cout << endl;

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
	//	allocate memory for particles
	// particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
	particle_t *particles = new particle_t[n];
	//	define a particle
	MPI_Datatype PARTICLE;
	MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
	MPI_Type_commit( &PARTICLE );
	//	set up the data partitioning across processors
	//	define the offsets
	int particle_per_proc = (n + n_proc - 1) / n_proc;
	// int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
	int* partition_offsets = new int[n_proc+1];
	for( int i = 0; i < n_proc+1; i++ )
		partition_offsets[i] = min( i * particle_per_proc, n );
	//	define the sizes
	// int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
	int* partition_sizes = new int[n_proc];
	for( int i = 0; i < n_proc; i++ )
		partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
	//	allocate storage for local partition
	int nlocal = partition_sizes[rank];
	particle_t* local = new particle_t[nlocal];
	//	initialize
	//	every process creates a pointer to a renderwindow
	sf::RenderWindow window;
	// create a circle object, must be visible to all proc
	sf::CircleShape circle(4, 8);
	circle.setFillColor( sf::Color::Black );
	if( rank == 0 )
	{
		srand( time( NULL ) );
		//	init the particles
		init_particles( n, particles );
		//	only one process creates a window
		// function call defined above
		create_window(window);
	}
	//	distribute the particles
	MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
	// color proc 0 particles blue
	for( int i = 0; i < nlocal; ++i)
	{
		local[i].color = sf::Color::Blue;
	}
	//	simulate a number of time steps
	double simulation_time = read_timer( );
	for( int step = 0; step < nSteps; step++ )
	// while( 1 )
	{
		//	collect all global data locally
		MPI_Allgatherv( local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );
		if( rank == 0 )
		{
			window.clear( sf::Color::White );
			for( int i = 0; i < n; ++i )
			{
				// scale the default radius by the mass
				// circle.setRadius(defaultRadius * massToDrawRadius * particles[i].mass);
				circle.setRadius(particles[i].displaySize);
				circle.setFillColor(particles[i].color);
				// fmod is floating point modulus
				// don't do modulus TODO
				circle.setPosition( fmod( particles[i].x * WINDOW_WIDTH, WINDOW_WIDTH), fmod(particles[i].y * WINDOW_HEIGHT, WINDOW_HEIGHT ) );
				// cout << particles[i].x << " " << particles[i].y << endl;
				// circle.setPosition( .5 * WINDOW_WIDTH, .5 * WINDOW_HEIGHT );
				window.draw( circle );
			}
			for( int i = 0; i < nlocal; ++i )
			{
				// scale the default radius by the mass
				// circle.setRadius(defaultRadius * massToDrawRadius * particles[i].mass);
				circle.setRadius(local[i].displaySize);
				circle.setFillColor(sf::Color::Red);
				// fmod is floating point modulus
				// don't do modulus TODO
				circle.setPosition( fmod( local[i].x * WINDOW_WIDTH, WINDOW_WIDTH), fmod(local[i].y * WINDOW_HEIGHT, WINDOW_HEIGHT ) );
				// cout << particles[i].x << " " << particles[i].y << endl;
				// circle.setPosition( .5 * WINDOW_WIDTH, .5 * WINDOW_HEIGHT );
				window.draw( circle );
			}
			/*
			sf::Image frame = window.capture();
			ostringstream stream;
			stream << "frame" << step << ".png";
			string frame_name = stream.str();
			frame.saveToFile(frame_name);
			*/
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
		{
			move( local[i] );
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

void init_particles( int n, particle_t* p )
{
	for( int i = 0; i < n; i++ ) 
	{
		p[i].x = 0.0;
		p[i].y = 0.0;
		p[i].vx = 0.0;
		p[i].vy = 0.0;
		p[i].ax = 0.0;
		p[i].ay = 0.0;
		p[i].mass = 1.0;
		p[i].displaySize = p[i].mass * 4.0;
		p[i].vx = rand() / (double)RAND_MAX * initialVelocity;
		p[i].vy = rand() / (double)RAND_MAX * initialVelocity;
		p[i].x = rand() / (double)RAND_MAX;
		p[i].y = rand() / (double)RAND_MAX;
		p[i].color = sf::Color::Black;
	}
}

void apply_force( particle_t &particle, particle_t &neighbor)
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

void move( particle_t &p )
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
