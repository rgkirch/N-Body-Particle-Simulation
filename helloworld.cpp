#include <stdio.h>
#include "omp.h"
int main() {
	// executes the next part in parallel
	#pragma omp parallel
	if( 1 )
		printf( "hello\n" );
	printf( "world\n" );
	//printf( "%d %d\n", omp_get_thread_num(), omp_get_num_threads() );
	// the threads will execute the next part one at a time
	//#pragma omp critical
	// garuentees mutual exclusion for changes in memory locations
	//#pragma omp atomic

}
