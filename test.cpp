#include <iostream>

int main()
{
int i=0;
int j=1;

int n = 5;
	for(i;i<(n-1);i++)
	{
		for(j=(i+1);j<n;j++)
		{
			std::cout << i << " " << j << "\n";
		}
	}


}
