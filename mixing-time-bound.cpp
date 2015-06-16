
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <string.h>
#include <time.h>
#include <cmath>

#include <algorithm>

#include <queue>
#include <vector>
#include <map>


using namespace std;


std::vector<double> calcBound(double epsilon, double n, double mi){

	std::vector<double> bounds;
	bounds.push_back((mi/(2-2*mi))*log(1/(2*epsilon)));
	bounds.push_back( (log(n) + log(1/epsilon))/(1-mi) );
    
    return bounds;
}

int main(int argc, char const *argv[])
{

	double n, mi;
	cout << "Enter n : " << endl;
	cin >> n;

	//cout << "Enter e : " << endl;
	//cin >> e;

	cout << "Enter mi : " << endl;
	cin >> mi;

	std::vector<double> bounds = calcBound(1/n, n, mi);
	cout << bounds[0] << " < T(epsilon) < " << bounds[1] << endl;

}