#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include "common-gravity.h"

using namespace std;

//	benchmarking program
int main( int argc, char **argv )
{
	int navg, nabsavg=0;
	double dmin, absmin=1.0,davg,absavg=0.0;
	double rdavg,rdmin;
	int rnavg;

	//	process command line parameters
	if( find_option( argc, argv, "-h" ) >= 0 )
	{
		printf( "Options:\n" );
		printf( "-h to see this help\n" );
		printf( "-n <int> to set the number of particles\n" );
		printf( "-o <filename> to specify the output file name\n" );
		printf( "-s <filename> to specify a summary file name\n" );
		printf( "-no turns off all correctness checks and particle output\n");
		return 0;
	}

	int n = read_int( argc, argv, "-n", 1000 );
	char *savename = read_string( argc, argv, "-o", NULL );
	char *sumname = read_string( argc, argv, "-s", NULL );

	//	set up MPI
	int n_proc, rank;
	MPI_Init( &argc, &argv );
	// ask the mpi communicator for n_proc and rank
	MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	// define the right and left neighbor
	int rightNeighbor = (rank + 1) % n_proc;
	int leftNeighbor = (rank + n_proc - 1) % n_proc;

	//	allocate generic resources
	// FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
	// FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;

	particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );

	// define a particle as 6 contiguous doubles
	MPI_Datatype PARTICLE;
	MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
	MPI_Type_commit( &PARTICLE );

	//	set up the data partitioning across processors
	int particle_per_proc = (n + n_proc - 1) / n_proc;
	int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
	for( int i = 0; i < n_proc+1; i++ )
		partition_offsets[i] = min( i * particle_per_proc, n );

	int max_partition_size = 0;

	int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
	for( int i = 0; i < n_proc; i++ )
	{
		partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
		if (partition_sizes[i] > max_partition_size)
		{
			max_partition_size = partition_sizes[i];
		}
	}


	//	allocate storage for local partition
	int nlocal = partition_sizes[rank];
	// particle_t *local = (particle_t*) malloc( nlocal * sizeof(particle_t) );
	particle_t *local = (particle_t*) malloc( nlocal* sizeof(particle_t) );
	particle_t *oneLocal = (particle_t*) malloc( max_partition_size* sizeof(particle_t) );
	particle_t *twoLocal = (particle_t*) malloc( max_partition_size* sizeof(particle_t) );

	//  initialize and distribute the particles (that's fine to leave it unoptimized)
	set_size( n );
	if( rank == 0 )
	{
		init_particles( n, particles );
	}
	MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );

	int oneLength = nlocal;
	int twoLength = 0;

	//	simulate a number of time steps
	double simulation_time = read_timer( );
	for( int step = 0; step < NSTEPS; step++ )
	// for( int step = 0; step < 1; step++ )
	{
		navg = 0;
		dmin = 1.0;
		davg = 0.0;

		// set all the local ones to be 0
		for(int i = 0; i < nlocal; ++i)
		{
			local[i].ax = local[i].ay = 0;
		}

		// we should initialize oneLocal with a copy of the local
		// oneLength = nlocal;
		// left = (particle_t*) malloc( oneLength* sizeof(particle_t) );
		// for (int k = 0; k < nlocal; ++k)
		// {
			// left[k] = local[k];
		// }

		// we will pass these around in a circle so we pass it as many times as we have procssors
		// each iteration passes it twice so we have iter < n_proc/2
		for (int iteration = 0; iteration < n_proc/2; ++iteration)
		{
			// MPI_Barrier(MPI_COMM_WORLD);
			// data, count, type, destination, tag, communicator
			// the even one begins the send and the odd one revieves
			if (rank % 2 == 0) {
				if (step == 0) {
					MPI_Send(&nlocal, 1, MPI_INT, rightNeighbor, 0, MPI_COMM_WORLD);
					MPI_Send(local, nlocal, PARTICLE, rightNeighbor, 0, MPI_COMM_WORLD);
				} else {
					// update local based on the data just passed in
					MPI_Send(&oneLength, 1, MPI_INT, rightNeighbor, 0, MPI_COMM_WORLD);
					MPI_Send(oneLocal, oneLength, PARTICLE, rightNeighbor, 0, MPI_COMM_WORLD);
				}
				// MPI_Barrier(MPI_COMM_WORLD);
				MPI_Recv(&oneLength, 1, MPI_INT, leftNeighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(oneLocal, oneLength, PARTICLE, leftNeighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				for (int iLocal = 0; iLocal < nlocal; ++iLocal)
					for (int iOne = 0; iOne < oneLength; ++iOne)
						apply_force( local[iLocal], oneLocal[iOne] );

				MPI_Send(&oneLength, 1, MPI_INT, rightNeighbor, 0, MPI_COMM_WORLD);
				MPI_Send(oneLocal, oneLength, PARTICLE, rightNeighbor, 0, MPI_COMM_WORLD);
				// MPI_Barrier(MPI_COMM_WORLD);
				MPI_Recv(&oneLength, 1, MPI_INT, leftNeighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(oneLocal, oneLength, PARTICLE, leftNeighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// MPI_Barrier(MPI_COMM_WORLD);
				for (int iLocal = 0; iLocal < nlocal; ++iLocal)
					for (int iOne = 0; iOne < oneLength; ++iOne)
						apply_force( local[iLocal], oneLocal[iOne] );
			// rank is odd
			} else {
				MPI_Recv(&oneLength, 1, MPI_INT, leftNeighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(oneLocal, oneLength, PARTICLE, leftNeighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// MPI_Barrier(MPI_COMM_WORLD);
				if (step == 0) {
					MPI_Send(&nlocal, 1, MPI_INT, rightNeighbor, 0, MPI_COMM_WORLD);
					MPI_Send(local, nlocal, PARTICLE, rightNeighbor, 0, MPI_COMM_WORLD);
				} else {
					MPI_Send(&twoLength, 1, MPI_INT, rightNeighbor, 0, MPI_COMM_WORLD);
					MPI_Send(twoLocal, twoLength, PARTICLE, rightNeighbor, 0, MPI_COMM_WORLD);
				}
				for (int iLocal = 0; iLocal < nlocal; ++iLocal)
					for (int iOne = 0; iOne < oneLength; ++iOne)
						apply_force( local[iLocal], oneLocal[iOne] );

				MPI_Recv(&twoLength, 1, MPI_INT, leftNeighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(twoLocal, twoLength, PARTICLE, leftNeighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// MPI_Barrier(MPI_COMM_WORLD);
				MPI_Send(&oneLength, 1, MPI_INT, rightNeighbor, 0, MPI_COMM_WORLD);
				MPI_Send(oneLocal, oneLength, PARTICLE, rightNeighbor, 0, MPI_COMM_WORLD);
				// MPI_Barrier(MPI_COMM_WORLD);
				for (int iLocal = 0; iLocal < nlocal; ++iLocal)
					for (int iTwo = 0; iTwo < twoLength; ++iTwo)
						apply_force( local[iLocal], twoLocal[iTwo] );
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);

		if( find_option( argc, argv, "-no" ) == -1 )
		{
			MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
			MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
			if (rank == 0){
				if (rnavg) {
					absavg +=  rdavg/rnavg;
					nabsavg++;
				}
				if (rdmin < absmin) absmin = rdmin;
			}
		}

		for( int i = 0; i < nlocal; i++ )
			move( local[i] );

		// MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );

		// TODO add a gather so particles is update to visualize

		// if (rank == 0) {
		//	   MPI_Gather(local, partition_sizes, PARTICLE, n, PARTICLE, 0, MPI_COMM_WORLD);

		//
		//	collect all global data locally (not good idea to do)
		//
		// MPI_Allgatherv( local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );

		//
		//	save current step if necessary (slightly different semantics than in other codes)
		//

		// if( find_option( argc, argv, "-no" ) == -1 )
		//	 if( fsave && (step%SAVEFREQ) == 0 )
		//	   save( fsave, n, particles );
	}
	simulation_time = read_timer( ) - simulation_time;

	if (rank == 0) {
	  printf( "n = %d, simulation time = %g seconds", n, simulation_time);

	  if( find_option( argc, argv, "-no" ) == -1 )
	  {
		if (nabsavg) absavg /= nabsavg;
	  //
	  //  -the minimum distance absmin between 2 particles during the run of the simulation
	  //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
	  //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
	  //
	  //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
	  //
	  printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
	  if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
	  if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
	  }
	  printf("\n");

	  // Printing summary data
	  // if( fsum)
		// fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
	}

	//	release resources
	// if ( fsum )
		// fclose( fsum );
	free( partition_offsets );
	free( partition_sizes );
	free( local );
	free( oneLocal );
	free( twoLocal );
	free( particles );
	// if( fsave )
		// fclose( fsave );

	MPI_Finalize( );

	return 0;
}
