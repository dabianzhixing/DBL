#include "DBL.h"
#include <iostream>
#include <boost/archive/binary_iarchive.hpp>
#include <time.h>
#include <sys/time.h>
#include <climits>
using namespace std;

int main(int argc, char **argv) {

	if (argc != 3) {
		cerr << "Argument Error\n" << endl;
		exit(EXIT_FAILURE);
	}

	printf("%s\n", argv[1]);
	printf("Reading index\n");
	DBL dbl;
	ifstream ifs(argv[1]);
	_CHECK(ifs);
	boost::archive::binary_iarchive ia(ifs);
	ia >> dbl;
	printf("Reading query\n");
	ifstream is(argv[2]);
	_CHECK(!is.bad());
	vector<int> u;
	vector<int> v;
	for (int m, n ; is >> m >> n ; )
	{
		u.push_back(m);
		v.push_back(n);
	}
	int sz = (int)u.size();

	bool *result = new bool[sz];
	fill_n(result, sz, false);

	printf("Start Query\n");
  	dbl.query(
  			u,
  			v,
  			sz,
  			result
		  );



}
