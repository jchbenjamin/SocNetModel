#include <iostream>
#include <vector>
#include "boolGraphMatrix.hpp"

boolGraphMatrix::boolGraphMatrix()
{
	n = 0;
	m.reserve(5000);
}
//create from encoding
boolGraphMatrix::boolGraphMatrix(int size, std::vector< bool > v)
{
	n = size;
	int i = 0;
	for(i;i<n;i++)
	{
		std::vector<bool> line;
		line.reserve(n);
		int j = i+1;
		for(j;j<n;j++)
		{	
			line.push_back(v.at(;
			v.pop_back();
			m[j][i] = m[i][j];
		}
	} 
}

std::vector< bool > boolGraphMatrix::getEncoding()
{
	int i = 0;
	std::vector< bool > enc;
	for(i;i<n;i++)
	{
		int j = i+1;
		for(j;j<n;j++)
		{
			enc.push_back( m[i][j] );
		}
	}
	
	return enc;
}

bool boolGraphMatrix::getVal(int i, int j)
{
	return m[i][j];
}

void boolGraphMatrix::addEdge(int i, int j)
{
	m[i][j] = true;
	m[j][i] = true;
}

void boolGraphMatrix::delEdge(int i, int j)
{
	m[i][j] = false;
	m[j][i] = false;
}

int boolGraphMatrix::getSize()
{
	return n;
}

std::vector<int> boolGraphMatrix::getNeighborsOfNode(int i)
{
	std::vector<int> r;
	int j = 0;
	for(j; j < n; j++)
	{
		if(m[i][j] == true)
		{
			r.push_back(j);
		}
	}
	return r;
}

void boolGraphMatrix::printM()
{
	int i=0;

	for(i;i<n;i++)
	{
		int j = 0;
		std::cout << "{";
		for(j;j<n;j++)
		{
			if(m[i][j] == true){
				std::cout << "t" << ",";
			}
			if(m[i][j] == false){
				std::cout << "f" << ",";
			}
		}
		std::cout << "}\n";
	}

}
