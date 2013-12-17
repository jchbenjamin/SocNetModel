

//Represents an economic agent in social networking model
class gameNode
{
	public:
		gameNode(int); //boring default constructor
		double getMargUtil(int, int);
		void addX(int);
		void subX(int);
		double putUtil(int, int);
		float getCost(int);
		void calcTotal();
		double getTotal();

	private:
		int n; //number of share types ('colors' in color graph analogy)
		std::vector< float > b; //coefficient for utility function u(x) = -bx + c  usually less than one
		std::vector< float > c; //constant for utility function usually positive integer
		std::vector< float > p; //price as discount for each share type  0<p<1
		std::vector< double > x; //current utility stack for each share type.
		std::vector< double > u; //current actual utility as function of x
		double total; //total utility
};
//utility functions
//u(x) = -(b/2)*x^2 + cx
//u'(x) = -bx+c

//Represents environment for the game
/*class gameEnv
{

	public:
		gameEnv(int,int);
		//~gameEnv();
		
	private:
		int numNodes;
		int numColors;
		boolGraphMatrix m;
		std::vector< int > nodeChoice;
		std::vector< gameNode > nodes;

};
*/



	
