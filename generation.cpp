
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
