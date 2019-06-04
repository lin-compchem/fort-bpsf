// Defines size_t, FILE, and fpos_t
#include <stdio.h>
// This is where size_t is defined
#include <stdlib.h>
#include <curl/curl.h>
#include <iostream>
#include <string>

using namespace std;
//
//  libcurl variables for error strings and returned data

static char errorBuffer[CURL_ERROR_SIZE];
static std::string buffer;


static int writer(char *data, size_t size, size_t nmemb,
                  std::string *writerData)
{
  if(writerData == NULL)
    return 0;

  writerData->append(data, size*nmemb);

  return size * nmemb;
}

std::string bufferToString(char* buffer, int bufflen)
{

}
static bool init(CURL *&conn, char *url)
{
  CURLcode code;

  conn = curl_easy_init();

  if(conn == NULL) {
    fprintf(stderr, "Failed to create CURL connection\n");
    exit(EXIT_FAILURE);
  }

  code = curl_easy_setopt(conn, CURLOPT_ERRORBUFFER, errorBuffer);
  if(code != CURLE_OK) {
    fprintf(stderr, "Failed to set error buffer [%d]\n", code);
    return false;
  }

  code = curl_easy_setopt(conn, CURLOPT_URL, url);
  if(code != CURLE_OK) {
    fprintf(stderr, "Failed to set URL [%s]\n", errorBuffer);
    return false;
  }

//  code = curl_easy_setopt(conn, CURLOPT_FOLLOWLOCATION, 1L);
//  if(code != CURLE_OK) {
//    fprintf(stderr, "Failed to set redirect option [%s]\n", errorBuffer);
//    return false;
//  }

  code = curl_easy_setopt(conn, CURLOPT_WRITEFUNCTION, writer);
  if(code != CURLE_OK) {
    fprintf(stderr, "Failed to set writer [%s]\n", errorBuffer);
    return false;
  }

  code = curl_easy_setopt(conn, CURLOPT_WRITEDATA, &buffer);
  if(code != CURLE_OK) {
    fprintf(stderr, "Failed to set write data [%s]\n", errorBuffer);
    return false;
  }

  return true;
}
//
//  libcurl connection initialization
//
int main(int arg, char **argv)
{
	CURL *curl;
	CURLcode result;
	char outs[] = "{\"instances\": [1.0, 2.0, 5.0]}";
    char url[] = "http://localhost:8501/v1/models/half_plus_two:predict";
	
	// Always start off with init
	curl = curl_easy_init();

	if(curl) {
		if(!init(curl, url)){
            fprintf(stderr, "Connection init failed\n");
			exit(EXIT_FAILURE);
		}
		// Set the url
		//result = curl_easy_setopt(curl, CURLOPT_URL, url);
		//if(result != CURLE_OK){

		//}
		// Set the predictions string
	    result = curl_easy_setopt(curl, CURLOPT_POSTFIELDS, outs);
		if(result != CURLE_OK){

		}
        // Set the writefunction, this is called whenever chunks of data
		// are generated
		//result = curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writer);
		//if(result != CURLE_OK) {
        //    cout << "uhoh couldn't set writefunction";
		//}
        //result = curl_easy_setopt(curl, CURLOPT_WRITEDATA, &buffer);
        //if(result != CURLE_OK) {
        //  fprintf(stderr, "Failed to set write data [%s]\n", errorBuffer);
        //  return false;
        //}
		result = curl_easy_perform(curl);
        ///* Check for errors */ 
		//// "\'{\"instances\": [1.0, 2.0, 5.0]}\'\n"
		//// "{\"instances\": [1.0, 2.0, 5.0]}"
        //if(result != CURLE_OK)
        //  fprintf(stderr, "curl_easy_perform() failed: %s\n",
        //          curl_easy_strerror(result));
 
        ///* always cleanup */ 
        //curl_easy_cleanup(curl);

	}
	curl_global_cleanup();

	cout << "Hello World\n";
	cout << outs;
	cout << '\n';
	printf("%s\n",buffer.c_str());
	cout << "\n!!!!!!!\n";
	return 0;
}
