/*
 * This contains a driver to check if tensorflow serving is running on a 
 * computer
 */
#include "tfserving_utils.hpp"
#include <iostream>

using namespace std;


int main(int argc, char **argv) {
    // Check to make sure they put thet model name when calling the program
    if (argc != 2) {
      cerr << "Error: There should only be one argument passed to the program"
          << endl;
      cerr << "    It should be the name of the tensorflow model. I.e.:" << endl;
      cerr << "    \"./demo1 [model_name]\"" << endl;
      exit(EXIT_FAILURE);
    }
    //TFServer a(8504, argv[1]);
    TFServer a(8504);
}
