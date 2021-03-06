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


//test fitness -- RUN TO SEE IF EVOLUTIONARY ALGORITHM CONVERGES ANSWER SHOULD BE NEAR 2000
/*
double gameEnv::fitness( std::shared_ptr< gene > g)
{
	double fit = 2000.0;
	int i =0;
	int j = 0;
	for(i;i<(n-1);i++)
	{
		for(j=i+1;j<n;j++)
		{
			if(g->edge(i,j))
			{
				fit -= 1.0;
			}
		}
	}
}
*/

//First Experiment
/*
double gameEnv::fitness( std::shared_ptr< gene > g)
{
	double fit = 0.0;
	int i = 0;
	int j = 0;
	int k = 0;
	for(i;i<(n-1);i++)
	{
		for(j=i+1;j<n;j++)
		{
			if(g->edge(i,j))
			{
				k=0;
				for(k;k<numU;k++)
				{
					fit += ( nodes[i]->uAt(k) * nodes[j]->uAt(k) );
				}
						
			}
		}
	}
	if(fit<0.0)
		return 0.0;
	else
		return fit;
}
*/
/*
//Second Experiment
double gameEnv::fitness( std::shared_ptr< gene > g)
{
	double fit = 0.0;
	double thisfit = 0.0;
	int connects = 0;
	int i = 0;
	int j = 0;
	int k = 0;
	for(i=0;i<n;i++)
	{
		connects = 0;
		thisfit = 0.0;
		for(j=0;j<n;j++)
		{
			if(g->edge(i,j))
			{
				connects++;
				k=0;
				for(k;k<numU;k++)
				{
					thisfit += ( nodes[i]->uAt(k) * nodes[j]->uAt(k) );
				}
						
			}
		}
		if(connects > 3)
		{
			thisfit = 0.0;
		}
		fit += thisfit;
	}
	if(fit<0.0)
		return 0.0;
	else
		return fit;
}
*/

/*
//Third Experiment
double gameEnv::fitness( std::shared_ptr< gene > g)
{
	double fit = 0.0;
	double thisfit = 0.0;
	double thisscore = 0.0;
	int connects = 0;
	int i = 0;
	int j = 0;
	int k = 0;
	for(i=0;i<n;i++)
	{
		connects = 0;
		thisfit = 0.0;
		for(j=0;j<n;j++)
		{
			if(g->edge(i,j))
			{
				connects++;
				k=0;
				for(k;k<numU;k++)
				{
					thisfit += ( 100 - pow(nodes[i]->uAt(k) - nodes[j]->uAt(k),2.0) );
				}
						
			}
		}
		if(connects > 5)
		{
			thisfit = 0.0;
		}
		fit += thisfit;
	}
	if(fit<0.0)
		return 0.0;
	else
		return fit;
}
*/
//Fourth Experiment  Gamma Functions
double gameEnv::fitness( std::shared_ptr< gene > g)
{
	double fit = 0.0;
	double thisfit = 0.0;
	double thisscore = 0.0;
	int connects = 0;
	int i = 0;
	int j = 0;
	int k = 0;
	for(i=0;i< (n-1);i++)
	{
		connects = 0;
		thisfit = 0.0;
		for(j=i+1;j<n;j++)
		{
			if(g->edge(i,j))
			{
				connects++;
				k=0;
				for(k;k<numU;k++)
				{
					thisfit += ( nodes[i]->uAt(k) * nodes[j]->uAt(k));
				}
						
			}
		}
		k=0;
		for(k;k<connects;k++)
		{
			thisfit = pow(thisfit,0.7);
		}
		fit += thisfit;
	}
	if(fit<0.0)
		return 0.0;
	else
		return fit;
}

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
