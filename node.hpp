#ifndef __NODE_HPP__
#define __NODE_HPP__

#include <vector>
#include <stdio.h>
#include <algorithm>
#include <memory>
#include <random>
#include <vector>
#include <cmath>
#include <iostream>

#include "gene.hpp"

class gameNode
{
	public:

		gameNode(std::vector< double >); //N is number of utils
		gameNode(int, int type, double, double);
		gameNode(int, int type, double, double, std::default_random_engine);
		//~gameNode();
		double uAt(int);

		void printNode();

	private:
		int numU;
		std::vector< double > u;
};

class gameEnv
{
	public:
		gameEnv(int N, int u, int type, double, double); //number of nodes //number of utility items
		gameEnv(int N, int u, int type, double uMean, double uSD, std::default_random_engine generator);
		double fitness( std::shared_ptr< gene > );
		void printEnv();
		int getNumU();
		double getU(int i, int j);
		void checkIdenticals();
		
	private:
		int n; //number of Nodes
		int e; //number of Edges
		int numU; //number of utilities
		std::vector< std::shared_ptr< gameNode > > nodes;
		
		std::vector< double > weights;
		void affinity( std::shared_ptr< gene >, int); //get affinity of each node to another if connection exists
		
};



//u(x) = - ((b/2)*x)^2 + cx
//u'(x) = -(b)x + c
//
#endif
