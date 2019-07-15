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
using namespace rapidjson;
using namespace std;

void JSON_Message::setupDocument(string &json_string) {
    if (verbose) {
        cout << "Created JSON_Message " << my_name << endl;
        cout << "Parsing JSON string: " << endl << json_string << endl;
    }
    if (document.Parse(json_string.c_str()).HasParseError()) {
        parseFail();
    }
    return;
}
void JSON_Message::initialize(string &json_string) {
    setupDocument(json_string);
    parseDocument();
}
//
//Constructors
//
JSON_Message::JSON_Message() {
    cout << "IGOTHERE" << endl;
}
//
// This constructor allows for initializing with string that contains JSON
// data 
JSON_Message::JSON_Message(string &json_string, bool verbose_flag) {
    verbose = verbose_flag;
   initialize(json_string);
}

void JSON_Message::parseDocument() {
    cout << "hello" << endl;
}
