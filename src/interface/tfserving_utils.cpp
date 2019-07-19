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
/*
 * Convert a string to a character array.
 * This is not being used right now.
 */
char* TFServer::string_to_char (string str) {
    char *out = new char[str.length() + 1];
    strcpy(out, str.c_str());
    return out; 
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
