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
#include "tfs_messages.hpp"
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
}
//
// This constructor allows for initializing with string that contains JSON
// data 
JSON_Message::JSON_Message(string &json_string, bool verbose_flag) {
    verbose = verbose_flag;
   initialize(json_string);
}
//
// Exit the program if the parse fails
void JSON_Message::parseFail() {
    cerr << "Failure to parse JSON message for " << my_name << endl;
    exit(EXIT_FAILURE);
}
// Holder for the parse document routine
int JSON_Message::parseDocument() {
}

ModelVersionStatus::ModelVersionStatus(string &json_string, bool verbose_flag) {
    verbose = verbose_flag;
    initialize(json_string); 
}
// Parse the class' document object for the flowing variables:
// Version
// State
// Status
// Error Code
// Error Message
int ModelVersionStatus::parseDocument(){
    if (verbose) { cout << "Entering parseDocument" << endl;}
   // Check that there is a model version status field in the response.
   // Otherwise, the server is probably not up..
    if (!document.HasMember(mmvs)) {
        cerr << "Error, could not find field: " << mmvs
           << "in JSON response." << endl;
        cerr << "Either the TensorflowServing response has changed or the "
          << "model server is not running." << endl;
        exit(EXIT_FAILURE); 
    }

    // Get the model version
    if (!document[mmvs][0].HasMember(mversion)) {
        cerr << "Error, could not find \"version\" field in " << mmvs 
            << " array." << endl;
        exit(EXIT_FAILURE);
    }
    if (!document[mmvs][0][mversion].IsString()) {
        cerr << "Error parsing version" << endl;
        exit(EXIT_FAILURE);
    }
    version = document[mmvs][0][mversion].GetString();
    // Get the model state
    if (!document[mmvs][0].HasMember(mstate)) {
        cerr << "Error, could not find \""<<mstate<<"\" field in " << mmvs 
            << " array." << endl;
        exit(EXIT_FAILURE);
    } 
    if (!document[mmvs][0][mstate].IsString()) {
        cerr << "Error parsing model state" << endl;
        exit(EXIT_FAILURE);
    }
    state = document[mmvs][0][mstate].GetString();
    // Get the model error code
    if (!document[mmvs][0][mstatus].HasMember(merror_code)) {
        cerr << "Error, could not find \""<< merror_code<<"\" field in " << mmvs 
            << " array." << endl;
        exit(EXIT_FAILURE);
    } 
    if (!document[mmvs][0][mstatus][merror_code].IsString()) {
        cerr << "Error parsing Error Code" << endl;
        exit(EXIT_FAILURE);
    }
    error_code = document[mmvs][0][mstatus][merror_code].GetString();
    // Get the model error message
    if (!document[mmvs][0][mstatus].HasMember(merror_message)) {
        cerr << "Error, could not find \""<< merror_message<<"\" field in " << mmvs 
            << " array." << endl;
        exit(EXIT_FAILURE);
    } 
    if (!document[mmvs][0][mstatus][merror_message].IsString()) {
        cerr << "Error parsing Error Message" << endl;
        exit(EXIT_FAILURE);
    }
    error_message = document[mmvs][0][mstatus][merror_message].GetString();
    
    if (verbose) {
        printStatus();
    }
    if (verbose) { cout << "Leaving parseDocument" << endl;}
}
void ModelVersionStatus::printStatus() {
    cout << "MODEL STATUS:" << endl;
    cout << "Version: " << version << endl;
    cout << "State: " << state << endl;
    cout << "Error code: " << error_code << endl;
    cout << "Error message: " << error_message << endl;

}
void ModelVersionStatus::checkStatus() {
    bool ok = true;
    if (error_code != good_code) {
        cerr << "Error, model error code should be \""<<good_code<<"\". It is " << error_code;
        ok = false;
    }
    if (state != good_state) {
        ok = false;
        cerr << "Error, model state should be \"AVAILABLE\". It is " << state;
    }
    if (!ok) {
        cerr << "Model Version Status check failed";
        exit(EXIT_FAILURE);
    }
    if (verbose) { cout << "Model Status is GOOD!" << endl; }
}


EnGradMessage::EnGradMessage(string &json_string, bool verbose_flag) {
    verbose = verbose_flag;
    message = json_string;
    initialize(json_string); 
}
//
// This really doesn't do much except make sure that there is a 
// field with 'outputs' in it
//
// The rest of the heavy lifitng is done when the user requests
// an energy or gradient term
int EnGradMessage::parseDocument() {
    // Assert that there is an output!
    if (!document.HasMember(mout)) {
        cerr << "Error, could not read outputs in response from TensorFlow "
             << "server" << endl;
        exit(EXIT_FAILURE);
    }
    return 0;
}
//
// Return the energy output
//
//double EnGradMessage::parseEnergy() {
//    if (document[mout][0][0].IsDouble()) {
//        return document[mout][0][0].GetDouble();
//    } 
//    cerr << "Error parsing energy from server response" << endl;
//    printMessage();
//    exit(EXIT_FAILURE);
//}
//
// This just prints the server output
void EnGradMessage::printMessage() {
    cerr << "Message from server:" << endl;
    cerr << message << endl;
}
//
// Here we parse the basis, really it is an arbitrary array parser but we are
// limited by having to select the right document object
double EnGradMessage::parseEnergy() {
    const char *name = mener.c_str();
    checkOutput(name); 
    if (!document[mout][name].IsArray()) {
        cerr << "Energy with name " << name << "is not JSON array" << endl;
    }
    return document[mout][name][0][0].GetDouble(); 
} 
//
// Here we're just looking to make sure the output field has a given object
void EnGradMessage::checkOutput(const char *name) {
    if (!document[mout].HasMember(name)) {
        cerr << "Error looking for an output in response from "
             << "TensorFlow server" << endl;
        cerr << "Expected object with name: {outputs:{" << name << "}}"
            << endl;
        exit(EXIT_FAILURE);
    }
    return;
} 
//
// Here we parse the basis, really it is an arbitrary array parser but we are
// limited by having to select the right document object
//void EnGradMessage::parseArray(const char* name, int x_size, int y_size) {
//    if (!document[mout].HasMember(name)) {
//        cerr << "Error, could not read gradient outputs in response from "
//             << "TensorFlow server" << endl;
//        cerr << "Expected object with name: {outputs:{" << name << "}}"
//            << endl;
//        exit(EXIT_FAILURE);
//    }
//    if (!document[mout][name].IsArray()) {
//        cerr << "Basis with name " << name << "is not JSON array" << endl;
//    }
//} 
