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
