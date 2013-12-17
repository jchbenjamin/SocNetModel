#include <stdlib.h>
#include "util.hpp"

float randomFloat(float a, float b) {
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = b - a;
    float r = random * diff;
    return a + r;
}

double myFit()
{
	double f = (double) randomFloat(0.0,400.0);
	return f;
}
