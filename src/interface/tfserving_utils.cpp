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
// RapidJSON Headers
#include <rapidjson/rapidjson.h>
#include <rapidjson/document.h>
#include <rapidjson/prettywriter.h>
// Program includes
#include "tfs_json_messages.h"

using namespace std;
using namespace rapidjson;


static char errorBuffer[CURL_ERROR_SIZE];

 
// Writer function modified from libcurl tutorial
//
// It appends the output from the request to the curl buffer.
// It will remove previous output (ithink) but keep other output?
int writer(char *data, size_t size, size_t nmemb, string *buffer) {
    int result = 0;
    if (buffer != NULL) {
         buffer->append(data, size * nmemb);
         result = size * nmemb;
         }
    return result;
}

// Send a http GET request to the curl server URI
// Then check the resulting message to make sure the server is running
// as expected.
bool check_tf_server(CURL *curl, string curl_buffer) {
    CURLcode code;
    code = curl_easy_setopt(curl, CURLOPT_HTTPGET, 1L);
    code = curl_easy_perform(curl);
    if(code != CURLE_OK) {
        fprintf(stderr, "Failed to get TF server status [%s]\n", errorBuffer);
    }
    ModelVersionStatus status(curl_buffer, true);
    return true; 
}


static bool init(CURL *&conn, const char *url, string curl_buffer)
{
    CURLcode code;
  
    conn = curl_easy_init();
  
    if(conn == NULL) {
      fprintf(stderr, "Failed to create CURL connection\n");
      exit(EXIT_FAILURE);
    }
    //Error buffer
    //code = curl_easy_setopt(conn, CURLOPT_ERRORBUFFER, errorBuffer);
    //if(code != CURLE_OK) {
    //  fprintf(stderr, "Failed to set error buffer [%d]\n", code);
    //  return false;
    //}
    //Server address
    code = curl_easy_setopt(conn, CURLOPT_URL, url);
    if(code != CURLE_OK) {
      fprintf(stderr, "Failed to set URL [%s]\n", errorBuffer);
      return false;
    }
    // Set the writefunction, this is called whenever chunks of data
    // are generated
    code = curl_easy_setopt(conn, CURLOPT_WRITEFUNCTION, writer);
    if(code != CURLE_OK) {
      fprintf(stderr, "Failed to set writer [%s]\n", errorBuffer);
      return false;
    }
    // This is where the output from the request goes.
    code = curl_easy_setopt(conn, CURLOPT_WRITEDATA, &curl_buffer);
    if(code != CURLE_OK) {
      fprintf(stderr, "Failed to set write data [%s]\n", errorBuffer);
      return false;
    }
    // Set the time out to a generous 2 seconds
    code = curl_easy_setopt(conn, CURLOPT_TIMEOUT, 2L);
    if(code != CURLE_OK) {
        fprintf(stderr, "Failed to set timeout option [%s]\n", errorBuffer);
    }
    return true;
}

string get_model_uri(int port, string model_name) {
    string uri = "http://localhost:";
    uri += to_string(port);
    uri += "/v1/models/" + model_name;
    return uri;
}

char* string_to_char (string str) {
    char *out = new char[str.length() + 1];
    strcpy(out, str.c_str());
    return out; 
}


int main(int arg, char **argv) {
    CURL *curl = NULL;
    CURLcode result;
    int port = 8500;
    bool inited = false;
    string model_name = "half_plus-two";
    string uri = get_model_uri(port, model_name) ;
    string curl_buffer; 
    cout << uri << "\n";

    // Initialize curl libraries
    curl_global_init(CURL_GLOBAL_DEFAULT);

    // Initialize connection with specific url
    inited = init(curl, uri.c_str(), curl_buffer);
    if (!inited) {
        cerr << "FAILED to INIT";
        exit(EXIT_FAILURE);
    }
    //
    // Check that the tensorflow server is OK
    check_tf_server(curl, curl_buffer);
    // Cleanup curl
    curl_global_cleanup();

}
