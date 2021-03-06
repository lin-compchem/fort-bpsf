/*
 * This file contains the function definitions for the capi for the TFServer
 * class
 */
#include "../parameters.h"
#include "tfserving_utils.hpp"
#include <iostream>


using namespace std;

//Constructor
TFSERVER* create_tfserver(int port, char *model_name) {
//TFSERVER* create_tfserver(int port) {
  #ifdef DEBUG
    cout << "C API, Create TFServer Object" << endl;
  #endif 
    return new TFServer(port, model_name);     
    //return new TFServer(port);     
}

//Destructor
void delete_tfserver(TFSERVER* tfserver) {
  #ifdef DEBUG
    cout << "C API, Delete TFServer Object" << endl;
  #endif 
    delete tfserver;
}
// Get the energy from the tensorflow server
void tfs_bpsf_energy(TFSERVER* tfserver,  double *basis, int *max_bas,
        int *max_atom, int *num_bas, int *num_atom, int *num_el,
        double *energy) {
  #ifdef DEBUG
    cout << "C API, Call sendBPSF Energy" << endl;
  #endif 
    return tfserver->sendBPSF(
            basis, max_bas, max_atom, num_bas, num_atom, num_el, energy
            );
}
// Get the energy from the tensorflow server
void tfs_bpsf_gradient(TFSERVER* tfserver,  double *basis, int *max_bas,
        int *max_atom, int *num_bas, int *num_atom, int *num_el,
        double *energy, double *gradient) {
  #ifdef DEBUG
    cout << "C API, Call sendBPSF Gradient" << endl;
  #endif 
    return tfserver->sendBPSF(
            basis, max_bas, max_atom, num_bas, num_atom, num_el, energy,
            gradient);
}
// This goes with 
void tfs_model_test1(TFSERVER* tfserver) {
    return tfserver->ModelTest1();
}
 
