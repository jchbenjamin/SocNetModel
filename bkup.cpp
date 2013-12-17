#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include "game.hpp"

gameNode::gameNode(int j)
{
	int i = 0;
	n = j;
	for(i;i<j;i++)
	{
		b.push_back(0.5);
		c.push_back(5.0);
		p.push_back(0.1/n);
		x.push_back(0.0);
		u.push_back(0.0);
	}
}

void gameNode::addX(int i)
{
	x[i] = std::min(x[i] + 1.0, (double) (c[i]/b[i]) );
}
void gameNode::subX(int i)
{
	x[i] = std::max(x[i] - 1.0, 0.0);
}
float gameNode::getCost(int i)
{
	if(i>=n) {
		std::cerr << "Error out of bounds gameNode::getCost";
		return 0.0;
	}
	return c[i];
}

double gameNode::getMargUtil(int i, int o) //o is offset
{
	if(i>=n) {
		std::cerr << "Error out of bounds gameNode::getMargUtil";
		return 0.0;
	}
	double m = -(double)(b[i]) * (x[i] + o) + (double)(c[i]); 
	return std::max(0.0, m);
}

double gameNode::putUtil(int i, int o)
{
	if(i>=n) {
		std::cerr << "Error out of bounds gameNode::putUtil";
		return 0.0;
	}
	
	double t = -((double)(b[i]/2)) * ( pow((double)(x[i] + o),2) ) + (double)(c[i]) * x[i];
	return t;
}

void gameNode::calcTotal()
{
	int i = 0;
	double tot= 0.0;
	for(i;i<n;i++)
	{
		tot += u[i];
	}
	total = tot;
}
double gameNode::getTotal()
{
	return total;
}
/*

gameEnv::gameEnv(int nodeInput, int colorInput)
{
	numNodes = nodeInput;
	numColors = colorInput;
	int j=0;
	for(j;j<numNodes;j++)
	{
		nodeChoice.push_back(0);
		gameNode g = gameNode(numColors);
		nodes.push_back(g);
	}
}



gameEnv::~gameEnv()
{
	int i = 0;
	for(i;i<numNodes;i++)
	{
		delete(*nodes[i]);
	}
}
*/	



