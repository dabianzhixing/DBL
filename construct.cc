#include "DBL.h"
#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <time.h>
#include <sys/time.h>

using namespace std;

int main(int argc, char **argv) {
  if (argc != 3) {
    printf("Argument Error");
    exit(EXIT_FAILURE);
  }

  printf("%s\n", argv[1]);
  DBL dbl;
  dbl.construct_index(argv[1]);

  printf("Serializing the Index......\n\n\n");
  ofstream ofs(argv[2],ios_base::binary);
  _CHECK(ofs);
  boost::archive::binary_oarchive oa(ofs);
  oa << dbl;
  ofs.close();

  exit(EXIT_SUCCESS);
}
