#include <iostream>
#include <map>

std::pair< int, int > nodesFromEdge(int t, int numEdges) 
{
	std::pair< int, int > ret;
	if(t >= numEdges)
	{
		std::cerr << "nodesFromEdge error, target must be less than total number of edges!\n";
		ret = std::make_pair(0,1);
		return ret;
	}
	int p = numNodes - 1;
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


int main()
{
	int numNodes = 20;
	int w = numNodes * (numNodes -1) / 2;
	int i = 0;
	for(i=0;i<w;i++)
	{
		std::pair<int,int> myPair = nodesFromEdge(i,numNodes);
		std::cout << "a: " << myPair.first << " b: " << myPair.second << "\n";
	}
}
