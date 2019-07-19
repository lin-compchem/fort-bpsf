/*
 * This file contains the function definitions for the capi for the TFServer
 * class
 */
#include "tfserving_utils.hpp"
#include <iostream>

#define DEBUG

using namespace std;

//Constructor
//TFSERVER* create_tfserver(int port, char *model_name) {
TFSERVER* create_tfserver(int port) {
  #ifdef DEBUG
    cout << "C API, Create TFServer Object" << endl;
  #endif 
    //return new TFServer(port, model_name);     
    return new TFServer(port);     
}

//Destructor
void delete_tfserver(TFSERVER* tfserver) {
  #ifdef DEBUG
    cout << "C API, Delete TFServer Object" << endl;
  #endif 
    delete tfserver;
}
