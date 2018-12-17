#include "main.hpp"

using namespace std;

int main(int argc, char **argv) {

  if (argc == 2) {
    char *infilename = argv[1];
    test(infilename, "out.root");
  } else if (argc == 3){
    char *infilename = argv[1];
    char *outfilename = argv[2];
    test(infilename, outfilename);

  }


  return 0;
}
