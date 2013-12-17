#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <new>
#include <random>
#include "util.hpp"
#include "gene.hpp"
#include "generation.hpp"
#include "node.hpp"

int main(int argc, char* argv[])
{
	srand(time(NULL));

/*
	int n = 10000;  //nodes in graph
	int p = 40;	//population size of pool
	int u = 10; //utilities
	int e = 5; //elites
	float r = 0.1; //mutation rate
	int c = 25; //crossovers
	int iterations = 7500;
*/

	std::vector< char > myCode = {1,0,1,1,1,0,0,1,0,1};
	gene geneOne = gene(5,myCode);
	geneOne.printCode();

//	geneOne.printCycles();


	return 0;
}
