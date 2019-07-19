/* 
 * This file contains "CURLER" classes which are objects that have their own
 * libcurl object, error buffer, etc.
 *
 * Note that CURL_GLOBAL_INIT must be called before calling these routines.
 */
#include <string>
#include <iostream>
#include "curler.hpp"
/*
 * Initialize the curl object.
 *
 * Here we are currently using the EASY version of CURL and are setting
 *  CURLOPT_ERRORBUFFER
 *  CURLOPT_WRITEFUNCTION
 *  CURLOPT_WRITER
 *  CURLOPT_TIMEOUT
 *  CURLOPT_URL
 */

using namespace std;

Curler::Curler() {
    // Main init 
    curl = curl_easy_init();
    if(curl == NULL) {
      cerr << "Failed to create CURL connection" << endl;
      exit(EXIT_FAILURE);
    }
    //Error buffer
    code = curl_easy_setopt(curl, CURLOPT_ERRORBUFFER, errorBuffer);
    if(code != CURLE_OK) {
      cerr << "Failed to set error buffer. Code: " <<  code << endl;
      exit(EXIT_FAILURE);
    }
    // Set the writefunction, this is called whenever chunks of data
    // are generated
    code = curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writer);
    if(code != CURLE_OK) {
      cerr << "Failed to write function. Error message: " <<  errorBuffer;
      exit(EXIT_FAILURE);
    }
    //// This is where the output from the request goes.
    code = curl_easy_setopt(curl, CURLOPT_WRITEDATA, &curl_buffer);
    if(code != CURLE_OK) {
      cerr << "Failed to set write data. Error message:" << errorBuffer
          << endl;
      exit(EXIT_FAILURE);
    }
    // Set the time out to a generous 2 seconds
    code = curl_easy_setopt(curl, CURLOPT_TIMEOUT, 2L);
    if(code != CURLE_OK) {
        cerr << "Failed to set timeout option. Error message: " << errorBuffer
            << endl;
      exit(EXIT_FAILURE);
    }
}
//
//Here we need to clean up the curl operation and stuff
//
Curler::~Curler() {
    curl_easy_cleanup(curl);
}
// Constructer that can also set URL!
//Curler::Curler(std::string uri){
//    Curler();
//    setURL(uri);
//}
/*
 * Set the URL for this curler
 */
void Curler::setURL(string url) {
    code = curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    if(code != CURLE_OK) {
      cerr << "Failed to set URL: " << url << endl << 
          " Error buffer:" << errorBuffer << endl;
      exit(EXIT_FAILURE);
    }
}
/* Writer function modified from libcurl tutorial
 *
 * It appends the output from the request to the curl buffer.
 * It will remove previous output (ithink) but keep other output?
 */
int Curler::writer(char *data, size_t size, size_t nmemb, string *buffer) {
    int result = 0;
    if (buffer != NULL) {
         buffer->append(data, size * nmemb);
         result = size * nmemb;
         }
    return result;
}

void Curler::httpGet() {
    code = curl_easy_setopt(curl, CURLOPT_HTTPGET, 1L);
    code = curl_easy_perform(curl);
    if(code != CURLE_OK) {
        cerr <<  "Failed to get TF server status \n Error buffer: " << errorBuffer;
        exit(EXIT_FAILURE);
    }
}
