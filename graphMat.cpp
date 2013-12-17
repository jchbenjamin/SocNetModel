#include <iostream>
#include <vector>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "boolGraphMatrix.hpp"

graphMat::graphMat(int size, std::vector< bool > v)
{
	n = size;
	m.resize(n-1, n-1);
	
	
	
	
