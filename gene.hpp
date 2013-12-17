#ifndef _GENE_HPP_
#define _GENE_HPP_

#include <iostream>
#include <vector>
#include <memory>
#include <stdlib.h>
#include <map>

#include "util.hpp"
//#include "gtfit.hpp"
//#include <boost/utility.hpp>

class vertex //structure for DFS
{
	public:
		int c; //color 0=white 1=gray 2=black
		int p; //pi
		int d; //discovery time
		int f; //finish time
};

class gene
{
	public:
		gene(int); //random creation
		gene(int, std::vector< char > );//create automatic code to test.
		gene( const std::shared_ptr< gene >, const std::shared_ptr< gene >, float); //crossover creation
		void mutate(float); //mutation randomly at constructor from parents
		void putFitness(double); //put fitness into gene, but check if >= 0.0
		double getFitness();
		std::vector< std::vector< int > > getCycles();
		void printCode();
		void printCycles(); //test to see if cycles algorithm works.
		void printFit();
		void printMask(std::vector< char >);
		std::vector< char > getCode();
		int getN();
		
		int numCycles(int i, int j);


	
		std::pair< int, int > nodesFromEdge(int, int);
		bool edge(int, int);
		friend class generation;


		
	private:
		int n; //NUMBER OF NODES
		std::vector< char > code; //representation of edges
		//std::vector< vertex > nodeList; //nodes structure for DFS
		std::vector< std::vector< char > > cycles;
		double fitness;
		
		
		void cycleExpand();
		std::vector< char > Xor (std::vector< char> a, std::vector< char > b);
		std::vector< char > convertToGene( std::vector< int > in); 
		bool edgeCycle(int a, int b, std::vector< char > in);
			
		void DFS();
		void DFSvisit(std::vector< vertex >* l, int* time, int u, std::vector< int >);
		
		
};



#endif // _GENE_HPP_
