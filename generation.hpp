#ifndef _GENERATION_HPP_
#define _GENERATION_HPP_

#include <vector>
#include <memory>
#include <fstream>
#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>

#include "node.hpp"
#include "gene.hpp"
#include "util.hpp"

class generation
{
	public:
		generation(int, int, int, int, float, std::shared_ptr< gameEnv > );
		generation(int, int, int, int, float, std::shared_ptr< gameEnv >, std::vector< std::shared_ptr< gene > > );
		std::vector< std::shared_ptr< gene > > select(); //pass gene pool to next generation //use in constructor for next generation
		
		void printPool();
		void printFits();
		void printGexf();
		void sort();
		
	private:
		int numNodes;
		int popSize;
		int elites;
		int numCrossOvers;
		float mutationRate;	
	
		void runFitness();
		
		std::shared_ptr< gameEnv > env;
		std::vector< std::shared_ptr< gene > > pool;
		std::shared_ptr< gene > chooseParent();
		
};

#endif
