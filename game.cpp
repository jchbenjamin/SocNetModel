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

	if(argc!=11)
	{
		std::cout << "nodes, population, utilities, elites , mut rate , crossovers, iterations, dist type(1 normal, 2 gamm), uMean, uStd\n";
		return 0;
	}
/*
	int n = 10000;  //nodes in graph
	int p = 40;	//population size of pool
	int u = 10; //utilities
	int e = 5; //elites
	float r = 0.1; //mutation rate
	int c = 25; //crossovers
	int iterations = 7500;
*/
	int n = atoi(argv[1]);
	int p = atoi(argv[2]);
	int u = atoi(argv[3]);
	int e = atoi(argv[4]);
	float r = atof(argv[5]);
	int c = atoi(argv[6]);
	int iterations= atoi(argv[7]);
	int dist = atoi(argv[8]);
	float mean = atof(argv[9]);
	float stddev = atof(argv[10]);
	
	

	std::random_device rd;
 	std::default_random_engine generator(rd());

	std::shared_ptr< gameEnv > gE = std::shared_ptr< gameEnv >(new gameEnv(n, u, dist, mean, stddev)); //criteria for environment //1 is normal dist //2 is gamma
	gE->printEnv();
	
	std::cout << "start!\n";
	std::cout << "size of gene vector:" << (n-1)*(n)/2 << "\n";
	std::cout << "genesize * population:" << ((n-1)*(n)/2)*p << "\n";

	generation g = generation(n,p,e,c,r,gE);
	//g.printPool();
	
	int i = 0;
	for(i;i< iterations ;i++)
	{
		g = generation(n,p,e,c,r,gE,g.select());
		std::cout << "generation " << i << "\n";
//		g.printPool();
	}
	
		std::cout << "Graph:: Nodes: " << n << " Utility rankings: " << u << "\n";
		std::cout << "Evolution:: Poolsize: " << p << " Elites: " << e << " Crossovers: " << c << " Mutation rate: " << r << "\n";
		g.sort();
		g.printFits();
		g.printGexf();
//	g.sort();
	//g.printFits();

	//g.printPool();
/*	
	gene geneOne = gene(10);
	geneOne.printCode();
	gene geneTwo = gene(10);
	geneTwo.printCode();
	
	gene geneMix = gene(geneOne.getCode(),geneTwo.getCode(),0.9);

	geneMix.printCode();
*/

	return 0;
}
