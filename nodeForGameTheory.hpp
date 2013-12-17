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

class gameNode
{
	public:

		gameNode(int N, std::vector< double >, std::vector< double >, std::vector< double >, std::vector< double >);
		gameNode(int, double, double, double, double, double, double, double, double);
		~gameNode();
		
		//void delUtil(int i, int xOffset);
		double margUtil(int i, double xOffset);

		void printNode();

		void turn(int i, std::vector< double > utils);

	private:
		int n;
		int nodeNum;

		void discount(int i);//apply cost of move to all x update totUtil;	
		void computeTotUtil(); //constructor helper function to initialize utility. After this, utility is incrementally changed.
		std::vector< double >* b;
		std::vector< double >* c;
		std::vector< double >* util;
		std::vector< double >* cost;

		double utilX(int i);

		double totUtil;

};

//u(x) = - ((b/2)*x)^2 + cx
//u'(x) = -(b)x + c
//
#endif
