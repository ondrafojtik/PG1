#include "stdafx.h"
#include "tutorials.h"

int main()
{
	printf( "PG1, (c)2011-2020 Tomas Fabian\n\n" );

	_MM_SET_FLUSH_ZERO_MODE( _MM_FLUSH_ZERO_ON );
	_MM_SET_DENORMALS_ZERO_MODE( _MM_DENORMALS_ZERO_ON );

	//return tutorial_1();
	//return tutorial_2();
#if LEGO_SPACESHIP
	return tutorial_3("../../../data/6887_allied_avenger.obj");
#endif

#if REFRACTION_SPHERE
	return tutorial_3("../../../data/geosphere.obj");
#endif

#if CORNELL_BOX
	return tutorial_3("../../../data/cornell_box2.obj");
#endif
}
