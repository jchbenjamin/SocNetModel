#include "node.hpp"

gameNode::gameNode(int N, std::vector< double > inB, std::vector< double > inC, std::vector< double > inUtil, std::vector< double > inCost)
{
	n = N;
	b = new std::vector< double >(inB);
	c = new std::vector< double >(inC);
	util = new std::vector< double >(inUtil);
	cost = new std::vector< double >(inCost);
}


void gameNode::computeTotUtil()
{
	totUtil=0.0;
	double temp=0.0;
	int i=0;
	for(i=0;i<n;i++)
	{
		temp = pow( -( (b->at(i)/2) * util->at(i) ), 2.0) + (c->at(i) * util->at(i));
		totUtil += temp;
	}
}

gameNode::gameNode(int N, double bMean, double bSD, double cMean, double cSD, double uMean, double uSD, double costMean, double costSD)
{

	n = N;
	int i = 0;
	std::vector< double > tempV;

	std::default_random_engine generator;
	std::normal_distribution<double> bDist(bMean,bSD);
	std::normal_distribution<double> cDist(cMean,cSD);
	std::normal_distribution<double> uDist(uMean,uSD);
	std::normal_distribution<double> costDist(costMean,costSD);

	tempV.push_back( 0.0 ); //non-participation b
	for(i=1;i<n;i++)
	{
		tempV.push_back( bDist(generator) );
	}
	b = new std::vector< double >(tempV);
	tempV.clear();

	tempV.push_back( 0.0 ); //non-participation c
	for(i=1;i<n;i++)
	{
		tempV.push_back( cDist(generator) );
	}
	c = new std::vector< double >(tempV);
	tempV.clear();

	tempV.push_back( 0.0 ); //utility of non-participation
	for(i=1;i<n;i++)
	{
		double tempU = uDist(generator);
		if(tempU > (c->at(i) / b->at(i)))
			tempU = (c->at(i) / b->at(i));
		else if(tempU < 0.0)
			tempU = 0.0;

		tempV.push_back( tempU );
	}
	util = new std::vector< double >(tempV);
	tempV.clear();

	tempV.push_back( 0.0 );
	for(i=1;i<n;i++)
	{
		double tempCost = costDist(generator);
		if(tempCost < 0.0)
			tempCost = 0.0;
		else if (tempCost > 1.0)
			tempCost = 1.0;

		tempV.push_back( tempCost );
	}
	cost = new std::vector< double >(tempV);

	computeTotUtil();
}


double gameNode::margUtil(int i, int xOffset)
{
	double maxU = c->at(i) / b->at(i);
	double marg = -(b->at(i))*(util->at(i) + (double)xOffset) + c->at(i);
	if(marg > maxU)
		marg = maxU;
	else if(marg < 0.0)
		marg = 0.0;

	return marg;
}

void gameNode::printNode()
{
	int i = 0;
	for(i=0;i<n;i++)
	{
		std::cout << "u" << i << ": " << util->at(i) << " ";
	}
		std::cout << "\n";
	for(i=0;i<n;i++)
	{
		std::cout << "b" << i << ": " << b->at(i) << " ";
	}
		std::cout << "\n";
	for(i=0;i<n;i++)
	{
		std::cout << "c" << i << ": " << c->at(i) << " ";
	}
		std::cout << "\n";

}
void gameNode::discount(int i)
{
	int j = 0;
	for(j;j<n;j++)
	{
		util->assign(j, util->at(j) * (1- cost->at(i)));
	}
}

void gameNode::turn(int i, std::vector< double > utils)
{

	int j = 0;
	for(j;j<n;j++)
	{
		util->assign(j, util->at(j) + utils[j]);
	}
	discount(i);

}

/*
void gameNode::delUtil(int i, int xOffset);
{
	tempUtil = -(( (b[utilNum]/2) * x[utilNum] ) ^2) + (c[utilNum] * x[utilNum]);
	totUtil -= tempUtil;

	if(util


}
*/


gameNode::~gameNode()
{
	delete b;
	delete c;
	delete util;
	delete cost;
}
