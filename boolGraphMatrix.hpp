#include <vector>

//represents undirected graph edges in matrix form
class boolGraphMatrix
{
	
	public:
		boolGraphMatrix(int, std::vector< bool>);
		boolGraphMatrix();

		//boolGraphMatrix(int,double,double, double) random variances of certain utilities
		bool getVal(int, int);
		void addEdge(int, int);
		void delEdge(int, int);
		int getSize();
		std::vector< bool > getEncoding();
		std::vector< int > getNeighborsOfNode(int);  //TODO PASS BY REFERENCE
		
		void printM();

		
	private:
		//outside is rows inside is columns for M[i][j] notation analogous to matrix notation conventions in math
		//simultaneous access to bool vector not thead safe but this data structure will only be read after construction
		std::vector< std::vector < bool > > m;
		//size of matrix and number of nodes
		int n;
};
