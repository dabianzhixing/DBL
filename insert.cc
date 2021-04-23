#include "DBL.h"
#include <fstream>
#include <boost/archive/binary_iarchive.hpp>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <iostream>

using namespace std;

int main(int argc, char **argv) {
	if (argc != 4)
	{
		cerr << "Argument Error" << endl;
		exit(EXIT_FAILURE);
	}

	printf("Reading index......\n");
	printf("%s\n",argv[1]);
	printf("Number of Insertions:%d\n", atoi(argv[3]));
	DBL dbl;
	ifstream ifs(argv[1]);
	_CHECK(ifs);
	boost::archive::binary_iarchive ia(ifs);
	ia >> dbl;

	dbl.insert_edges(argv[2], atoi(argv[3]));
}
