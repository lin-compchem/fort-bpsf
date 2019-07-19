/* 
 * This file contains "CURLER" classes which are objects that have their own
 * libcurl object, error buffer, etc.
 *
 * Note that CURL_GLOBAL_INIT must be called before calling these routines.
 */
#ifndef CURLER
#define CURLER
#include <curl/curl.h>
#include <string>
class Curler {
  public:
    // Variables
    char errorBuffer[CURL_ERROR_SIZE] = {' '};
    std::string curl_buffer; 
    CURL *curl = NULL;
    CURLcode code; // Return code for curl operations
    // Methods
    Curler();
    ~Curler();
    // Curler(std::string url); Cant get to work...
    void setURL(std::string url); 
    void httpGet();
  private:
    static int writer(char *data, size_t size, size_t nmemb, std::string *buffer);
    //bool init(); 
};
#endif 
