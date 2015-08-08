#include <iostream>

#include "constants.h"
#include "gem.h"
#include "matrix.h"
#include "vector.h"
#include "quaternion.h"


void main()
{
	auto m = gem::matrix4f{ 1, 1, 1, 1,
						    2, 2, 2, 2,
						    3, 3, 3, 3,
						    4, 4, 4, 4 };
	
	std::cout << "m is \n" << m << std::endl << "trace is " << gem::trace(m);
	
	getchar();
}