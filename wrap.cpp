#include <stdio.h>

inline float wrap( float num, float mod )
{
	return num - 2.0 * mod * (float)(int)( num / 2.0 * mod );
}

int main() {
	float five = 5.0;
	printf( "wrap: %g\n", wrap( -5.0, 4.0 ) );
}

