/* This file contains utility functions for working with tensorflow-serving
 * models.
 *
 * It can query a server and check the status
 */
#include <curl/curl.h>
#include <string>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
// Program includes
#include "tfs_messages.hpp"
#include "tfserving_utils.hpp"
#include "curler.hpp"

using namespace std;
using namespace rapidjson;
// Constructor
//TFServer::TFServer(int set_port, char *set_model_name)
//                       :port(set_port), model_name(set_model_name) {
TFServer::TFServer(int set_port, char* my_model_name):port(set_port) {
//TFServer::TFServer(int set_port):port(set_port) {
   model_name = my_model_name;
   start_interface();  
}
// Destructor
TFServer::~TFServer() {
    end_interface();
}
//
// This returns the URL for the model status, has information about whether
// the model is running and if there is an error code.
string TFServer::get_model_status_uri() {
    string uri = "http://localhost:";
    uri += to_string(port);
    uri += "/v1/models/" + model_name;
    return uri;
}
//
// This returns the URI for the metadata, which has information about
// what the name of the inputs and outputs for a given model are.
string TFServer::get_metadata_uri() {
    string uri = "http://localhost:";
    uri += to_string(port);
    uri += "/v1/models/" + model_name + "/metadata";
    return uri;
}
/* Send a http GET request to the curl server URI
 * Then check the resulting message to make sure the server is running
 * as expected.
 */
void TFServer::check_tfs_status() {
    // Local variables
    string uri = get_model_status_uri();
    Curler stat_curl;
    
    // Initialize connection with the model url
    stat_curl.setURL(uri.c_str());
    stat_curl.httpGet();

    // Check the return from the curl operation 
    ModelVersionStatus status(stat_curl.curl_buffer, true);
    status.checkStatus();
}
/* Send a http GET request to the tensorflow serving metadata URL
 */
void TFServer::get_tfs_metadata() {
    // Local variables
    string uri = get_metadata_uri();
    Curler meta_curl;
    
    // Initialize connection with the model url
    meta_curl.setURL(uri.c_str());
    meta_curl.httpGet();

    // Check the return from the curl operation 
    cout << meta_curl.curl_buffer; 
    //ModelVersionStatus status(meta_curl.curl_buffer, true);
    //status.checkStatus();
}
void TFServer::sendBPSF(double *basis, int *max_bas, int *max_atom,
                   int *num_bas, int *num_atom, int *num_el, double *energy) {
    cout << "SENT BASIS" << endl;
    print_basis(basis, max_bas, max_atom, num_bas, num_atom, num_el);
    exit(EXIT_FAILURE);
}

void TFServer::start_interface() {
    // Initialize curl libraries
    curl_global_init(CURL_GLOBAL_DEFAULT);
    
    // Check that the tensorflow server is OK
    check_tfs_status();
    
    // Get the model metadata
    get_tfs_metadata();
} 

void TFServer::end_interface() {
    // Cleanup curl
    curl_global_cleanup();
}
// Prints the basis
// There are a lot of parameters here
// The basis object has [Fortran style] indices of
//  (max_bas, max_atom, total_els). However, many of these will be blank.
// Therefore we need to only index over the number of things that we have.
// This is done by the counters num_bas, numatom, and num_el
//
// double *basis: array of size (max_bas, max_atom, num_el)
//                this is the basis to print
void TFServer::print_basis(double *basis, int *max_bas, int *max_atom,
        int *num_bas, int *num_atom, int *num_el){
    int i, j;
    for (i=0; i < *num_el; ++i)
    {
        j = i * *max_bas * *max_atom;
        print_el_basis(&basis[j], max_bas[0], max_atom[0], num_bas[i],
               num_atom[i]);
    }
}
void TFServer::print_el_basis(double *basis, int max_bas, int max_atom,
        int num_bas, int num_atom) {
    int i, j, idx;
    int basis_counter = 0;
    for (i=0; i < num_atom; ++i) {
        basis_counter = i * max_bas;
        for (j=0; j < num_bas; ++j){
            idx = basis_counter + j;
            printf("%f           %d \n", basis[idx], idx);
        }
    }
}


