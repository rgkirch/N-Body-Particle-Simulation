#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdarg.h>    
#include <vector>
#include <sys/time.h>
#include <SDL.h>
#include <SDL_opengl.h>
#include <GL/glu.h>

#define WINDOW_SIZE 1
#define DEFAULT_FILENAME "data.txt"

struct particle_t
{
	float x, y;
};
	
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized ) {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

void make_surface()
{
	// draw points with proper filtering, otherwise draw aliased points
    glEnable( GL_POINT_SMOOTH );
	// glend the computed fragment color values with the values in the color buffers
    glEnable( GL_BLEND );
	// something about blending the colors coming in with the colors already there
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glPointSize( 2 );
    glClearColor( 1, 1, 1, 1 );	
    glViewport (0, 0, WINDOW_SIZE, WINDOW_SIZE);
	// specifies which matrix stack is the target for subsequent matrix operations
	// projection - apply subsequent operations to the matrix stack
    glMatrixMode(GL_PROJECTION);
	// replaces the current matrix with the identity matrix
    glLoadIdentity();  
}

int main( int argc, char *argv[] )
{
    const char *filename = argc > 1 ? argv[1] : DEFAULT_FILENAME;
    FILE *f = fopen( filename, "r" );
    if( f == NULL )
    {
        printf( "failed to find %s\n", filename );
        return 1;
    }
    int n;
	int count = 0;
    fscanf( f, "%d", &n );

	double* particles = (double*) malloc( n * 2 * sizeof(double));
	double tempPair[2] = {0.0, 0.0};
	
	// if it is false that you have any more values to read in the current frame
	// AND
	// if it is also false that there is another pair of values to begin another frame
	// THEN
	// you are done so exit
	// OTHERWISE
	// there are more values to read in the current frame
	while (1)
	{
		while(fscanf( f, "%g%g", &p.x, &p.y ) == 2)
		{
			particles.push_back( p );
    fclose( f );
	
    // int nframes = particles.size( ) / n;
    // if( nframes == 0 )
        // return 2;
    
    // int window_size = (int)((size+2*eps)*SCALE);
    // window_size = window_size > MIN_SIZE ? window_size : MIN_SIZE;
	
    // init graphics system
    if( SDL_Init( SDL_INIT_VIDEO ) != 0 ) {
        printf( "SDL_Init failed: %s\n", SDL_GetError( ) );
        return 1;
    }
	SDL_Window *window =  SDL_CreateWindow( "Particle Simulation", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WINDOW_SIZE, WINDOW_SIZE, SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE );
    if( window == NULL ) {
        printf( "SDL_SetVideoMode failed: %s\n", SDL_GetError( ) );
        exit( 1 );
    }
    SDL_GLContext glcontext = SDL_GL_CreateContext(window);
    make_surface( );
    for( bool done = false; !done; )
    {
        SDL_Event event;
        while( SDL_PollEvent( &event ) )
        {
            if( ( event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_ESCAPE ) || ( event.type == SDL_QUIT ) )
                done = true;
            else if( event.window.event == SDL_WINDOWEVENT_RESIZED )
                make_surface( event.window.data1, event.window.data2 );
        }
		
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
		
        glMatrixMode( GL_PROJECTION );
        glLoadIdentity( );
        gluOrtho2D( 0, WINDOW_SIZE, 0, WINDOW_SIZE );
		
        glColor3f( 0.75, 0.75, 0.75 );
        glBegin( GL_LINE_LOOP );
        glVertex2d( 0, 0 );
        glVertex2d( size, 0 );
        glVertex2d( size, size );
        glVertex2d( 0, size );
        glEnd( );
		
        int iframe = (int)(read_timer()*FPS) % nframes;
        particle_t *p = &particles[iframe*n];
		
        glColor3f( 0, 0, 0 );
		// defines how it should treat the following data, up to glEnd()
		// GL_POINTS - treats each verted as a single point
        glBegin( GL_POINTS );
        for( int i = 0; i < n; i++ )
			// glVertex2fv() specifies a vertex
			// accepts a pointer to an array of two elemnts
			// p is a struct but the elements are contiguous regardless
            glVertex2fv( &p[i].x );
        glEnd( );
		
        SDL_GL_SwapWindow(window);
    }
	SDL_GL_DeleteContext(glcontext);
    SDL_Quit();
	
    return 0;
}
