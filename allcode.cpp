###compile.sh
#!/bin/bash

g++ game.cpp -c -g -std=c++0x -rdynamic
g++ gene.cpp -c -g -std=c++0x -rdynamic
g++ generation.cpp -c -g -std=c++0x -rdynamic
g++ node.cpp -c -g -std=c++0x -rdynamic
g++ util.cpp -c -g -std=c++0x -rdynamic

g++ -o game game.o gene.o generation.o node.o util.o -rdynamic -std=c++0x



/*********************************************
*  GAME.CPP
*********************************************/


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

/*********************************************
*  GENE.HPP
*********************************************/
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

/*********************************************
*  GENE.CPP
*********************************************/
#include "gene.hpp"
#include "util.hpp"

gene::gene(int nodes) //generate gene of matrix size n
{
	//	std::cout << "gene constructor: random!\n";
	n = nodes;
	
	int m = (n-1) * (n) / 2;
		code.reserve(m);
	int i = 0;
	for(i;i<m;i++)
	{
	
		if( (rand() % 2) == 0) //random number either zero or one represents random bit
		{
			code.push_back(1);
		} else
		{
			code.push_back(0);
		}
	}
	
	//PUT IN FITNESS FUNCTION
	fitness = 0.0;
	//DFS();
	//printCycles();
	//cycleExpand();
}

gene::gene(int inN, std::vector< char > inG)
{
	n = inN;
	
	fitness = 0.0;
	code = inG;
//	DFS();
//	printCycles();
	//cycleExpand();
}

gene::gene(const std::shared_ptr< gene > a, const std::shared_ptr< gene > b, float p)
{
		//std::cout << "gene constructor: go\n";
	if(a->code.size() != b->code.size())
	{
		std::cout << a->code.size() << "a size\n";
		std::cout << b->code.size() << "b size\n";
		std::cerr << "gene constructor: vector size mismatch\n";
		return;
	}
	n = a->getN();
	
	int m = (n-1) * n / 2;
	//std::cout << "size" << n << "\n";
	
	int cutOne = 0;
	int cutTwo = 0;
	while(cutOne == cutTwo)
	{
		cutOne = rand() % m;
	//	std::cout << "cutOne" << cutOne << "\n";
		cutTwo = rand() % m;
	//	std::cout << "cutTwo" << cutTwo << "\n";
	}
	if(cutOne <= cutTwo)
	{
		int i = 0;
		for (i;i<cutOne;i++)
		{
			code.push_back(a->code.at(i));
		}
		for (i;i<cutTwo;i++)
		{
			code.push_back(b->code.at(i));
		}
		for (i;i<n;i++)
		{
			code.push_back(a->code.at(i));
		}
	}else
	{
		int i = 0;
		for (i;i<cutTwo;i++)
		{
			code.push_back(b->code.at(i));
		}
		for (i;i<cutOne;i++)
		{
			code.push_back(a->code.at(i));
		}
		for (i;i<n;i++)
		{
			code.push_back(b->code.at(i));
		}
	}
	//random mutation of offspring
	mutate(p);
	
	//INIT FITNESS;
	fitness = 0.0;
//	DFS();
	//cycleExpand();
}


void gene::mutate(float p) //USE RANDOM MUTATE TECHNIQUE //p is percentage mutation
{

	int m = (n-1) * n / 2;

	std::vector< char > newCode;
	if(p < 0.0 || p > 1.0)
	{
		std::cerr << "incorrect mutation percentage";
		return;
	}
	std::vector< char > mask;
	int i = 0;
	for(i;i<m;i++) //generate mask
	{
		if(randomFloat(0.0,1.0) < p)
		{
			mask.push_back(1);
		} else
		{
			mask.push_back(0);
		}
	}
	//printMask(mask);
	i=0;
	for(i;i<m;i++) //apply mask
	{
		newCode.push_back( (code[i] + mask[i]) % 2 );
	}
	code.clear();
	code = std::vector< char >(newCode);
	return;
}
	
std::vector< char > gene::getCode()
{
	return code;
}
int gene::getN()
{
	return n;
}


void gene::putFitness(double in)
{
	if(in > 0.0)
		fitness = in;
	else
		fitness = 0.0;
}


double gene::getFitness()
{
	return fitness;
}


std::pair< int, int > gene::nodesFromEdge(int t, int numEdges) 
{
	std::pair< int, int > ret;
	if(t >= numEdges)
	{
		std::cerr << "nodesFromEdge error, target must be less than total number of edges!\n";
		ret = std::make_pair(0,1);
		return ret;
	}
	int p = numEdges- 1;
	int i = 0;
	int j = 1;
	int q = 1;
	while( (p-q) < t)
	{
		i +=1;
		j +=1;
		t -= (p-q+1);
		q++;
	}
	j+=t;	
	ret = std::make_pair(i,j);
}


bool gene::edge(int a, int b)
{
//	std::cout << "EDGE\n";
	int m = n-1;
	int c = 0;
	if(a == b)
	{
		return false;
	}#include "gene.hpp"
#include "util.hpp"

gene::gene(int nodes) //generate gene of matrix size n
{
	//	std::cout << "gene constructor: random!\n";
	n = nodes;
	
	int m = (n-1) * (n) / 2;
		code.reserve(m);
	int i = 0;
	for(i;i<m;i++)
	{
	
		if( (rand() % 2) == 0) //random number either zero or one represents random bit
		{
			code.push_back(1);
		} else
		{
			code.push_back(0);
		}
	}
	
	//PUT IN FITNESS FUNCTION
	fitness = 0.0;
	//DFS();
	//printCycles();
	//cycleExpand();
}

gene::gene(int inN, std::vector< char > inG)
{
	n = inN;
	
	fitness = 0.0;
	code = inG;
//	DFS();
//	printCycles();
	//cycleExpand();
}

gene::gene(const std::shared_ptr< gene > a, const std::shared_ptr< gene > b, float p)
{
		//std::cout << "gene constructor: go\n";
	if(a->code.size() != b->code.size())
	{
		std::cout << a->code.size() << "a size\n";
		std::cout << b->code.size() << "b size\n";
		std::cerr << "gene constructor: vector size mismatch\n";
		return;
	}
	n = a->getN();
	
	int m = (n-1) * n / 2;
	//std::cout << "size" << n << "\n";
	
	int cutOne = 0;
	int cutTwo = 0;
	while(cutOne == cutTwo)
	{
		cutOne = rand() % m;
	//	std::cout << "cutOne" << cutOne << "\n";
		cutTwo = rand() % m;
	//	std::cout << "cutTwo" << cutTwo << "\n";
	}
	if(cutOne <= cutTwo)
	{
		int i = 0;
		for (i;i<cutOne;i++)
		{
			code.push_back(a->code.at(i));
		}
		for (i;i<cutTwo;i++)
		{
			code.push_back(b->code.at(i));
		}
		for (i;i<n;i++)
		{
			code.push_back(a->code.at(i));
		}
	}else
	{
		int i = 0;
		for (i;i<cutTwo;i++)
		{
			code.push_back(b->code.at(i));
		}
		for (i;i<cutOne;i++)
		{
			code.push_back(a->code.at(i));
		}
		for (i;i<n;i++)
		{
			code.push_back(b->code.at(i));
		}
	}
	//random mutation of offspring
	mutate(p);
	
	//INIT FITNESS;
	fitness = 0.0;
//	DFS();
	//cycleExpand();
}


void gene::mutate(float p) //USE RANDOM MUTATE TECHNIQUE //p is percentage mutation
{

	int m = (n-1) * n / 2;

	std::vector< char > newCode;
	if(p < 0.0 || p > 1.0)
	{
		std::cerr << "incorrect mutation percentage";
		return;
	}
	std::vector< char > mask;
	int i = 0;
	for(i;i<m;i++) //generate mask
	{
		if(randomFloat(0.0,1.0) < p)
		{
			mask.push_back(1);
		} else
		{
			mask.push_back(0);
		}
	}
	//printMask(mask);
	i=0;
	for(i;i<m;i++) //apply mask
	{
		newCode.push_back( (code[i] + mask[i]) % 2 );
	}
	code.clear();
	code = std::vector< char >(newCode);
	return;
}
	
std::vector< char > gene::getCode()
{
	return code;
}
int gene::getN()
{
	return n;
}


void gene::putFitness(double in)
{
	if(in > 0.0)
		fitness = in;
	else
		fitness = 0.0;
}


double gene::getFitness()
{
	return fitness;
}


std::pair< int, int > gene::nodesFromEdge(int t, int numEdges) 
{
	std::pair< int, int > ret;
	if(t >= numEdges)
	{
		std::cerr << "nodesFromEdge error, target must be less than total number of edges!\n";
		ret = std::make_pair(0,1);
		return ret;
	}
	int p = numEdges- 1;
	int i = 0;
	int j = 1;
	int q = 1;
	while( (p-q) < t)
	{
		i +=1;
		j +=1;
		t -= (p-q+1);
		q++;
	}
	j+=t;	
	ret = std::make_pair(i,j);
}


bool gene::edge(int a, int b)
{
//	std::cout << "EDGE\n";
	int m = n-1;
	int c = 0;
	if(a == b)
	{
		return false;
	}
	else if (a > b)
	{
		int temp = a;
		a = b;
		b = temp;
	}
	else if (b > m)
	{
		std::cerr << "invalid edge call higher node out of bounds\n";
		return false;
	}

	for(int i=0;i<a;i++)
	{
		c+= m;
		m--;
	}
	c += ((b - a) - 1);
	if(code.at(c) == (char) 0)
	{
		return false;
	} else
	{
		return true;
	}
	
}

void gene::printFit()
{
	std::cout << fitness;
}

/*
void gene::printCycles()
{
	int i = 0;
	int j = 0;
	for(i;i<cycles.size();i++)
	{
	
		std::cout << "cycle: ";
		for(j=0;j<cycles[i].size();j++)
		{
			std::cout << cycles[i][j] << "->";
		}
		std::cout << "\n";
		
	}
}
*/


void gene::printCode()
{
	int i = 0;
	std::cout << fitness;
	std::cout << "{";
	for(i;i<code.size();i++)
	{
			if(code[i] == 1)
			{
				std::cout << "t,";
			} else
			{
				std::cout << "f,";
			}
	}
		std::cout << "}\n";
	
}

/*
void gene::printCycles()
{
	int i = 0;
	int j = 0;
	for(i;i<cycles.size();i++)
	{
		j=0;
	std::cout << i;
	std::cout << "{";
	for(j;j<cycles[i].size();j++)
	{
			if(cycles[i][j] == 1)
			{
				std::cout << "t,";
			} else
			{
				std::cout << "f,";
			}
	
	}
	std::cout << "}\n";
	}
}
*/

void gene::printMask(std::vector< char > mask)
{
	int i = 0;
	std::cout << "MSK";
	std::cout << "{";
	for(i;i<mask.size();i++)
	{
			if(mask[i] == 1)
			{
				std::cout << "t,";
			} else
			{
				std::cout << "f,";
			}
	}
		std::cout << "}\n";
	
}

void gene::DFS()
{

	std::vector< vertex >* l = new std::vector< vertex >();
	int u = 0;
	for(u;u<n;u++)
	{
		vertex thisV;
		thisV.c = 0;
		thisV.p = 0;
		thisV.d = 0;
		thisV.f = 0;
		l->push_back(thisV);
	}
	
	int time = 0;
	for(u=0;u<n;u++)
	{
//	std::cout << "dfsstart\n";
		if(l->at(u).c == 0)
		{
			//std::cout << "dfs" << u << "\n";
			std::vector< int > tree;
			tree.push_back(u);
			DFSvisit(l,&time, u, tree);
		}
	}
	
	delete l;
}


void gene::DFSvisit(std::vector< vertex >* l, int* time, int u, std::vector< int > branch)
{
	*time += 1;
	l->at(u).d = *time;
	l->at(u).c = 1; //gray
	int v = 0;
	for(v;v<n;v++)
	{
		//std::cout << "dfsvisit" << u << v << "\n";
		if(edge(u,v))
		{
			//std::cout << "edge" << u << v << "\n";
			if(l->at(v).c == 0)
			{
			//	std::cout << "testing\n";
				l->at(v).p = u;
				branch.push_back(v);
				DFSvisit(l,time,v,branch);
			}
			else if(l->at(v).c == 1 && branch.at(branch.size()-2) != v ) //we found a back edge so return path v->u in branch.
			{
				//std::cout << "cycle found!\n";
				int i = 0;
				while(branch[i] != v)
				{
					i++;
				}
				//if(branch[i] == v)
				//	std::cout << "success\n";
				
				std::vector< int > myCycle;
				for(i;i<branch.size();i++)
				{
			//		std::cout << "push " << branch[i] << "\n";
					myCycle.push_back(branch[i]);
				}
				std::cout<< "converting...";
				cycles.push_back(convertToGene(myCycle));
			//	std::cout << " I just walked away";
			}
		}
		//std::cout << "dfsvisitend" << u << v << "\n";
	}
	l->at(u).c = 2; //2 is black
	*time += 1;
	l->at(u).f = *time;
}

std::vector< char > gene::convertToGene( std::vector< int > in) {

	std::vector< char > outGene;
	


	int i = 0;
	
	for(i;i < (n*(n-1)/2); i++ )
	{
		outGene.push_back(0);
	}
	i=0;
	int a = 0;
	int b = 0;
	int m = n-1;
	int c = 0;
	for(i;i<( in.size()-1 );i++)
	{
		a = in.at(i);
		b = in.at(i+1);

//EDGE CODE
	c=0;

	if (a > b)
	{
		int temp = a;
		a = b;
		b = temp;
	}
	for(int j=0;j<a;j++)
	{
		c+= m;
		m--;
	}
	c += ((b - a) - 1);
	std::cout << "cycle bit at " << c << "\n";
	outGene.at(c) = 1;	
	/////EDGE CODE

	return outGene;
	
	}

}

std::vector< char > gene::Xor (std::vector< char> a, std::vector< char > b)
{
	std::vector< char > myCode;
	int i = 0;
	for(i;i < a.size();i++)
	{
		if(a.at(i) == 0 && b.at(i) == 0)
		{
			myCode.push_back(1);
		}
		else if(a.at(i) == 1 && b.at(i) == 1)
		{
			myCode.push_back(1);
		}
		else
		{
			myCode.push_back(0);
		}
	}
	return myCode;
}

void gene::cycleExpand()
{
	int i = 0;
	int j = 0;
	int s = cycles.size();
	for(i;i<(s-1);i++)
	{
		for(j = i+1;j<s;j++)
		{
			cycles.push_back(Xor(cycles.at(i) , cycles.at(j) ));
		}
	}

}

int gene::numCycles(int i, int j)
{
	int numba = 0;
	int p = cycles.size();
	int q = 0;
	for(q;q<p;q++)
	{
		if(edgeCycle(i,j,cycles[q]))
		{
			numba++;
		}
	}
	return numba;
}

bool gene::edgeCycle(int a, int b, std::vector< char > in)
{
//	std::cout << "EDGE\n";
	int m = n-1;
	int c = 0;
	if(a == b)
	{
		return false;
	}
	else if (a > b)
	{
		int temp = a;
		a = b;
		b = temp;
	}
	else if (b > m)
	{
		std::cerr << "invalid edge call higher node out of bounds\n";
		return false;
	}

	for(int i=0;i<a;i++)
	{
		c+= m;
		m--;
	}
	c += ((b - a) - 1);
	if(in.at(c) == (char) 0)
	{
		return false;
	} else
	{
		return true;
	}
	
}
	else if (a > b)
	{
		int temp = a;
		a = b;
		b = temp;
	}
	else if (b > m)
	{
		std::cerr << "invalid edge call higher node out of bounds\n";
		return false;
	}

	for(int i=0;i<a;i++)
	{
		c+= m;
		m--;
	}
	c += ((b - a) - 1);
	if(code.at(c) == (char) 0)
	{
		return false;
	} else
	{
		return true;
	}
	
}

void gene::printFit()
{
	std::cout << fitness;
}

/*
void gene::printCycles()
{
	int i = 0;
	int j = 0;
	for(i;i<cycles.size();i++)
	{
	
		std::cout << "cycle: ";
		for(j=0;j<cycles[i].size();j++)
		{
			std::cout << cycles[i][j] << "->";
		}
		std::cout << "\n";
		
	}
}
*/


void gene::printCode()
{
	int i = 0;
	std::cout << fitness;
	std::cout << "{";
	for(i;i<code.size();i++)
	{
			if(code[i] == 1)
			{
				std::cout << "t,";
			} else
			{
				std::cout << "f,";
			}
	}
		std::cout << "}\n";
	
}

/*
void gene::printCycles()
{
	int i = 0;
	int j = 0;
	for(i;i<cycles.size();i++)
	{
		j=0;
	std::cout << i;
	std::cout << "{";
	for(j;j<cycles[i].size();j++)
	{
			if(cycles[i][j] == 1)
			{
				std::cout << "t,";
			} else
			{
				std::cout << "f,";
			}
	
	}
	std::cout << "}\n";
	}
}
*/

void gene::printMask(std::vector< char > mask)
{
	int i = 0;
	std::cout << "MSK";
	std::cout << "{";
	for(i;i<mask.size();i++)
	{
			if(mask[i] == 1)
			{
				std::cout << "t,";
			} else
			{
				std::cout << "f,";
			}
	}
		std::cout << "}\n";
	
}

void gene::DFS()
{

	std::vector< vertex >* l = new std::vector< vertex >();
	int u = 0;
	for(u;u<n;u++)
	{
		vertex thisV;
		thisV.c = 0;
		thisV.p = 0;
		thisV.d = 0;
		thisV.f = 0;
		l->push_back(thisV);
	}
	
	int time = 0;
	for(u=0;u<n;u++)
	{
//	std::cout << "dfsstart\n";
		if(l->at(u).c == 0)
		{
			//std::cout << "dfs" << u << "\n";
			std::vector< int > tree;
			tree.push_back(u);
			DFSvisit(l,&time, u, tree);
		}
	}
	
	delete l;
}


void gene::DFSvisit(std::vector< vertex >* l, int* time, int u, std::vector< int > branch)
{
	*time += 1;
	l->at(u).d = *time;
	l->at(u).c = 1; //gray
	int v = 0;
	for(v;v<n;v++)
	{
		//std::cout << "dfsvisit" << u << v << "\n";
		if(edge(u,v))
		{
			//std::cout << "edge" << u << v << "\n";
			if(l->at(v).c == 0)
			{
			//	std::cout << "testing\n";
				l->at(v).p = u;
				branch.push_back(v);
				DFSvisit(l,time,v,branch);
			}
			else if(l->at(v).c == 1 && branch.at(branch.size()-2) != v ) //we found a back edge so return path v->u in branch.
			{
				//std::cout << "cycle found!\n";
				int i = 0;
				while(branch[i] != v)
				{
					i++;
				}
				//if(branch[i] == v)
				//	std::cout << "success\n";
				
				std::vector< int > myCycle;
				for(i;i<branch.size();i++)
				{
			//		std::cout << "push " << branch[i] << "\n";
					myCycle.push_back(branch[i]);
				}
				std::cout<< "converting...";
				cycles.push_back(convertToGene(myCycle));
			//	std::cout << " I just walked away";
			}
		}
		//std::cout << "dfsvisitend" << u << v << "\n";
	}
	l->at(u).c = 2; //2 is black
	*time += 1;
	l->at(u).f = *time;
}

std::vector< char > gene::convertToGene( std::vector< int > in) {

	std::vector< char > outGene;
	


	int i = 0;
	
	for(i;i < (n*(n-1)/2); i++ )
	{
		outGene.push_back(0);
	}
	i=0;
	int a = 0;
	int b = 0;
	int m = n-1;
	int c = 0;
	for(i;i<( in.size()-1 );i++)
	{
		a = in.at(i);
		b = in.at(i+1);

//EDGE CODE
	c=0;

	if (a > b)
	{
		int temp = a;
		a = b;
		b = temp;
	}
	for(int j=0;j<a;j++)
	{
		c+= m;
		m--;
	}
	c += ((b - a) - 1);
	std::cout << "cycle bit at " << c << "\n";
	outGene.at(c) = 1;	
	/////EDGE CODE

	return outGene;
	
	}

}

std::vector< char > gene::Xor (std::vector< char> a, std::vector< char > b)
{
	std::vector< char > myCode;
	int i = 0;
	for(i;i < a.size();i++)
	{
		if(a.at(i) == 0 && b.at(i) == 0)
		{
			myCode.push_back(1);
		}
		else if(a.at(i) == 1 && b.at(i) == 1)
		{
			myCode.push_back(1);
		}
		else
		{
			myCode.push_back(0);
		}
	}
	return myCode;
}

void gene::cycleExpand()
{
	int i = 0;
	int j = 0;
	int s = cycles.size();
	for(i;i<(s-1);i++)
	{
		for(j = i+1;j<s;j++)
		{
			cycles.push_back(Xor(cycles.at(i) , cycles.at(j) ));
		}
	}

}

int gene::numCycles(int i, int j)
{
	int numba = 0;
	int p = cycles.size();
	int q = 0;
	for(q;q<p;q++)
	{
		if(edgeCycle(i,j,cycles[q]))
		{
			numba++;
		}
	}
	return numba;
}

bool gene::edgeCycle(int a, int b, std::vector< char > in)
{
//	std::cout << "EDGE\n";
	int m = n-1;
	int c = 0;
	if(a == b)
	{
		return false;
	}
	else if (a > b)
	{
		int temp = a;
		a = b;
		b = temp;
	}
	else if (b > m)
	{
		std::cerr << "invalid edge call higher node out of bounds\n";
		return false;
	}

	for(int i=0;i<a;i++)
	{
		c+= m;
		m--;
	}
	c += ((b - a) - 1);
	if(in.at(c) == (char) 0)
	{
		return false;
	} else
	{
		return true;
	}
	
}

/*********************************************
*  GENERATION.HPP
*********************************************/
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

/*********************************************
*  GENERATION.CPP
*********************************************/

#include "gene.hpp"
#include "generation.hpp"
#include "util.hpp"

generation::generation(int m, int p, int e, int cross, float a, std::shared_ptr< gameEnv > envIn)
{

	numNodes = m;
	popSize = p;
	elites = e;
	numCrossOvers = cross;
	mutationRate = a;
	pool.reserve(popSize);
	int i = 0;
	
	//THROW EXCEPTIONS
	/*
	if(elites > matSize)
	  THROW
	 
	else if (numCrossOvers > (popSize-elites))
		THROW

	else if (0.0>mutationRate || 1.0 < mutationRate)
		THROW
	*/
	
	for(i;i<popSize;i++)
	{
		pool.push_back( std::shared_ptr< gene >(new gene(numNodes)) );
	}
	env = envIn;
	
	runFitness();
}

generation::generation(int m, int p, int e, int cross, float a, std::shared_ptr< gameEnv > envIn, std::vector< std::shared_ptr< gene > > g)
{
	numNodes = m;
	popSize = p;
	elites = e;
	numCrossOvers = cross;
	mutationRate = a;
	int i = 0;
	for(i=0;i<g.size();i++)
	{
		pool.push_back(g.at(i));
	}
	for(i;i<popSize;i++)
	{
		pool.push_back( std::shared_ptr< gene >(new gene(numNodes)) );
	}
	env = std::shared_ptr< gameEnv >(envIn);
	
	runFitness();
}
/*
bool generation::sortByFitness (std::shared_ptr< gene > i, std::shared_ptr< gene > j)
	{ return (i->getFitness() > j->getFitness() ); }
*/
std::shared_ptr< gene > generation::chooseParent()
{

	float myRoll = randomFloat(0.0,1.0); //rollin on roulette wheel
	std::vector< float > percentages;
	float percOffset = 0.0;
	int t = 0;  //total fitness
	int i = 0;
	for(i=0;i<popSize;i++)
	{
		t += pool[i]->getFitness();
	}
	
	for(i=0;i<popSize;i++)  //get roulette wheel percentages
	{
		percOffset += (float)pool[i]->getFitness() / t; 

		if(percOffset<1.0)
			{ percentages.push_back(percOffset); }
		else
			{ percentages.push_back(1.0); }
	}

	for(i=0;i<popSize;i++)
	{
		if(percentages[i] >= myRoll)
		{
			std::shared_ptr< gene > myChoice = pool[i];
			//pool.erase(pool.begin()+i);
			return myChoice;
		}
	}
}

/*generation::generation(int m, int p, float a, float b)
*/

bool sortByFitness (std::shared_ptr< gene > i, std::shared_ptr< gene > j)
	{ return (i->getFitness() > j->getFitness() ); }

std::vector< std::shared_ptr< gene > > generation::select()
{
		//std::cout << "select\n";
	std::vector< std::shared_ptr< gene > > nextPool;
	std::vector< std::shared_ptr< gene > > tempPool = std::vector< std::shared_ptr< gene > >(pool);

	//sort genes by fitness descending
	std::sort(tempPool.begin(),tempPool.end(),sortByFitness);
	int i = 0;

	//save elites

	for(i=0;i<elites;i++)
	{
		nextPool.push_back(tempPool[i]);
	}
	
	//next randomly breed new genes which may or may not mutate
	for(i=0;i<numCrossOvers;i++)
	{
	
		std::shared_ptr< gene > parentOne = chooseParent();
//		std::cout << "PARENT1   \n";
//		parentOne->printCode();
		
		std::shared_ptr< gene > parentTwo = chooseParent();
//		std::cout << "PARENT2   \n";
//		parentTwo->printCode();
			
		std::shared_ptr< gene > child(new gene( parentOne, parentTwo, mutationRate));
		
//		std::cout << "MYCROSSOVER   \n";
//		child->printCode();
//		std::cout << "push child\n";
		nextPool.push_back(child);
	}
		
	//delete tempPool;dreamy
	return nextPool;
	
}

void generation::sort()
{
	std::sort(pool.begin(),pool.end(),sortByFitness);
}
void generation::printFits()
{
	int i=0;
	for(i;i<popSize;i++)
	{
		pool[i]->printFit();
		std::cout << ", ";
	}
	std::cout << "\n";
}
void generation::printPool()
{
	int i = 0;
	for(i;i<popSize;i++)
	{
		pool[i]->printCode();
	}
	std::cout << "-----------------------\n";
}

void generation::runFitness()
{
	int i = 0;
	for(i;i<popSize;i++)
	{
//		std::cout << "runfit" << i << "\n";
		pool[i]->putFitness( env->fitness(pool.at(i)) );
	}		
}

void generation::printGexf()
{

	std::ofstream f;
	int u = env->getNumU();
	int i = 0;
	int j = 0;
	
	sort();
	std::shared_ptr< gene > g = std::shared_ptr< gene >(pool[0]);
	
	f.open("out.gexf");
	f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	f << "<gexf xmlns=\"http://www.gexf.net/1.1draft\" version=\"1.1\">\n";
   f << "<graph mode=\"static\" defaultedgetype=\"undirected\">\n";
   f << "<attributes class=\"node\">\n";

	for(i;i < u;i++)
	{
		f << "\t<attribute id=\"" << i << "\" title=\"u" << i << "\" type=\"float\">\n";
		f << "\t\t<default>0.0</default>\n";
		f << "\t</attribute>\n";
	}
	f << "</attributes>\n";
	f << "<nodes>\n";
	
	for(i=0;i<numNodes;i++)
	{
		f << "\t<node id=\"" << i << "\" label=\"" << i << "\" >\n";
		f << "\t<attvalues>\n";
		for(j=0;j<u;j++)
		{
			f << "\t\t<attvalue for=\"" << j << "\" value=\"" << env->getU(i,j) << "\" />\n";
		}
		f << "\t</attvalues>\n";
		f << "\t</node>\n";
	}
	f << "</nodes>\n";	
	f << "<edges>\n";
	for(i=0;i < (numNodes - 1) ;i++)  //logic from fitness function
	{
		for(j=(i+1);j < numNodes ;j++)
		{
			if(g->edge(i,j) == 1)
			{
				f << "\t<edge id=\"i" << i << "j" << j << "\" source=\"" << i << "\" target=\"" << j << "\" />\n";
			}
		}
	}
	f << "</edges>\n";
	f << "</graph>\n";
	f << "</gexf>\n";
	  
	f.close();
}


/*********************************************
*  NODE.HPP
*********************************************/
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


#endif

/*********************************************
*  NODE.CPP
*********************************************/
#include "node.hpp"

gameEnv::gameEnv(int N, int u, int type, double uMean, double uSD)
{
	n = N;
	numU = u;
	int i = 0;

    std::random_device rd;
  std::default_random_engine generator(rd());

	for(i;i<n;i++)
	{
		nodes.push_back( std::shared_ptr< gameNode >(new gameNode(u, type, uMean, uSD)) );
	}
	checkIdenticals();
}

gameEnv::gameEnv(int N, int u, int type, double uMean, double uSD, std::default_random_engine generator)
{
	n = N;
	numU = u;
	int i = 0;

	for(i;i<n;i++)
	{
		nodes.push_back( std::shared_ptr< gameNode >(new gameNode(u, type, uMean, uSD, generator)) );
	}
	checkIdenticals();
}

void gameEnv::checkIdenticals()
{
	int i = 0;
	int j = 1;
	for(i;i<(n-1);i++)
	{
		for(j;j<n;j++)
		{
			int uToTest = 0;
			int errCount = 0;
			for(uToTest; uToTest < numU; uToTest++)
			{
				if(getU(i,uToTest) == getU(j,uToTest))
				{
					errCount++;
				}
				if(errCount == (numU-1))
				{
					//ERROR WE FOUND IDENTICAL NODES
	std::cout << "ERROR!!!! IDENTICAL NODES!\n";
	std::cout << "ERROR!!!! IDENTICAL NODES!\n";
	std::cout << "ERROR!!!! IDENTICAL NODES!\n";
	std::cout << "ERROR!!!! IDENTICAL NODES!\n";
	std::cout << "ERROR!!!! IDENTICAL NODES!\n";
	std::cout << "ERROR!!!! IDENTICAL NODES!\n";
	std::cout << "ERROR!!!! IDENTICAL NODES!\n";
	std::cout << "ERROR!!!! IDENTICAL NODES!\n";
	std::cout << "ERROR!!!! IDENTICAL NODES!\n";			
				
				}
			}	
		}
	}
	return;

	
}
void gameEnv::printEnv()
{
	int i = 0;
	for(i;i<n;i++)
	{
		nodes[i]->printNode();
	}
}

/* THIRD FITNESS FUNCTION*/
/*
double gameEnv::fitness(std::shared_ptr< gene > g)
{
	affinity(g, 3); //affinity and max connects 3
	double fit = 0.0;
	int i = 0;
	int j=0;
	int q =0;
	int connects = 0;
	double thisfit=0.0;
	for(i;i < weights.size() ;i++) //iterate thru edges
	{
		fit += weights[i]; //weights generated by affinity function
	}
	
	i=0;
	
	for(i=0;i<n;i++)
	{
		q=0;
		thisfit = 0.0;
		connects = 0;
		for(j=0;j<n;j++)
		{
//			std::cout << "at " << i << " " << j;
			if(g->edge(i,j) == 1)
			{
				connects++;
				//std::cout << "connectiterate: " << connects << "\n";
				int k = 0;
				for(k;k<numU;k++)
				{
					thisfit = thisfit + ( nodes[i]->uAt(k) * nodes[j]->uAt(k) );
//										std::cout << "fit!" << fit << " " << k << "";
				}
				k=0;
//				std::cout << "add fit " << fit;
				for(k; k < g->numCycles(i,j); k++)
				{
					thisfit = pow(thisfit,2.0);
				}
			}
//			std::cout << "\n";
		}

	//	std::cout << "thisfit " << thisfit << "\n";
	//	std::cout << "conectmod " << (double)( n - connects) / (double) n << "\n";
		for(q; q<connects;q++)
		{
			thisfit = pow(thisfit,0.5);
		}

		
	//	std::cout << "corrected " << thisfit << "\n";
		fit = fit + thisfit;
	//	std::cout << "fit " << fit << "\n";	
	}	
	
	
	if(fit < 0.0)
	{
		fit = 0.0;
	}
	return fit;

}
*/

/* SECOND FITNESS FUNCTION*/

double gameEnv::fitness(std::shared_ptr< gene > g)
{
//	std::cout << "fitness func " << n << "\n";
	double fit = 0.0;  //n number of nodes
	int i = 0;
	int j = 0;
	int q = 0;
	int connects = 0;
	double thisfit = 0.0;
//	int m = n-1;
	for(i=0;i<n;i++)
	{
		q=0;
		thisfit = 0.0;
		connects = 0;
		for(j=0;j<n;j++)
		{
//			std::cout << "at " << i << " " << j;
			if(g->edge(i,j) == 1)
			{
				connects++;
				//std::cout << "connectiterate: " << connects << "\n";
				int k = 0;
				for(k;k<numU;k++)
				{
					thisfit = thisfit + ( nodes[i]->uAt(k) * nodes[j]->uAt(k) );
//										std::cout << "fit!" << fit << " " << k << "";
				}
//				std::cout << "add fit " << fit;
			}
//			std::cout << "\n";
		}

	//	std::cout << "thisfit " << thisfit << "\n";
	//	std::cout << "conectmod " << (double)( n - connects) / (double) n << "\n";
		for(q; q<connects;q++)
		{
			thisfit = pow(thisfit,0.5);
		}
	//	std::cout << "corrected " << thisfit << "\n";
		fit = fit + thisfit;
	//	std::cout << "fit " << fit << "\n";	
	}
	//std::cout << "endfit" << fit << "\n";
	return fit;
}


/* FIRST FITNESS FUNCTION */

/*
double gameEnv::fitness(std::shared_ptr< gene > g)
{
//	std::cout << "fitness func " << n << "\n";
	double fit = 0.0;  //n number of nodes
	int i = 0;
	int j = 1;
//	int m = n-1;
	for(i;i<(n-1);i++)
	{
		for(j=(i+1);j<n;j++)
		{
//			std::cout << "at " << i << " " << j;
			if(g->edge(i,j) == 1)
			{
				int k = 0;
				for(k;k<numU;k++)
				{
					fit = fit + ( nodes[i]->uAt(k) * nodes[j]->uAt(k) );
//										std::cout << "fit!" << fit << " " << k << "";
				}
//				std::cout << "add fit " << fit;
			}
//			std::cout << "\n";
		}	
	}
	//std::cout << "endfit" << fit << "\n";
	return fit;
}
*/

gameNode::gameNode(std::vector< double > inU)
{
	numU = inU.size();
	u = inU;
}

gameNode::gameNode(int N, int type, double uMean, double uSD)
{

	numU = N;
	int i = 0;
	std::vector< double > tempV;
   std::random_device rd;
 	std::default_random_engine generator(rd());

if(type==1)
{
//	std::gamma_distribution< double > uDist(uMean,uSD);
	std::normal_distribution<double> uDist(uMean,uSD);

	for(i=0;i<numU;i++)
	{
		tempV.push_back( uDist(generator) );
	}
	u = std::vector< double >(tempV);
} else if(type==2)
{
	std::gamma_distribution< double > uDist(uMean,uSD);
//	std::normal_distribution<double> uDist(uMean,uSD);

	for(i=0;i<numU;i++)
	{
		tempV.push_back( uDist(generator) );
	}
	u = std::vector< double >(tempV);
}

}

gameNode::gameNode(int N, int type, double uMean, double uSD, std::default_random_engine generator)
{

	numU = N;
	int i = 0;
	std::vector< double > tempV;

	std::gamma_distribution< double > uDist(uMean,uSD);
//	std::normal_distribution< double > uDist(uMean,uSD);

	for(i=0;i<numU;i++)
	{
		tempV.push_back( uDist( generator) );
	}
	u = std::vector< double >(tempV);
}

void gameNode::printNode()
{
	int i = 0;
	for(i=0;i<numU;i++)
	{
		std::cout << "u" << i << ": " << u.at(i) << " ";
	}
		std::cout << "\n";
}

double gameNode::uAt(int i)
{
	return u.at(i);
}

int gameEnv::getNumU()
{
	return numU;
}

double gameEnv::getU(int i, int j)
{
	return nodes[i]->uAt(j);

}
/*
void gameNode::delUtil(int i, int xOffset);
{
	tempUtil = -(( (b[utilNum]/2) * x[utilNum] ) ^2) + (c[utilNum] * x[utilNum]);
	totUtil -= tempUtil;

	if(util


}
*/

/*
gameNode::~gameNode()
{
	delete u;
}
*/


void gameEnv::affinity(std::shared_ptr< gene > g, int maxConnect)
{
	weights.clear(); //clear weights vector for new weights!
	int connectCount = 0;
	double affinity = 0.0;
	int k = 0;
	int i = 0;
	int j = 1;
//	int m = n-1;
	for(i;i<(n-1);i++)
	{
		connectCount = 0;
		for(j=(i+1);j<n;j++)
		{
			affinity = 0.0;
//			std::cout << "at " << i << " " << j;
			if(g->edge(i,j))
			{
				connectCount++;
				if(connectCount > maxConnect)
				{
					affinity = -10.0;
				} else
				{
					k = 0;
					for(k;k<numU;k++)
					{
						affinity += pow( nodes[i]->uAt(k) + nodes[j]->uAt(k), 2.0);
					}
//					std::cout << "add fit " << fit;
				}
			}
//			std::cout << "\n";
			weights.push_back(affinity);

		}

	}

}

/*********************************************
*  UTIL.HPP
*********************************************/
#ifndef _UTIL_HPP_
#define _UTIL_HPP_

#include <memory>
#include "gene.hpp"

float randomFloat(float, float);
double myFit();
#endif
/*********************************************
*  UTIL.CPP
*********************************************/
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

