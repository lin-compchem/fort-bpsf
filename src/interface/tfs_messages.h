#pragma once

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

using namespace rapidjson;
using namespace std;

#ifndef TFS_JSON
#define TFS_JSON
class JSON_Message {
    public:
        // General Vars
        Document document; // This is the object to parse the JSON string into
        const char my_name[8] = "generic"; // Name of message for error messages
        // Member functions
        JSON_Message();
        JSON_Message(string &json_string, bool verbose_flag=false);
    protected:
        bool verbose = false; // Flag for verbosity
        void setupDocument(string &json_string); 
        void parseDocument();
        void parseFail() {
            cerr << "Failure to Parse " << my_name << " Message" << endl;
            exit(EXIT_FAILURE);
        }
        void initialize(string &json_string);
};
#endif
#ifndef TFS_MVS
#define TFS_MVS
class ModelVersionStatus {
    public:
        int a;
        ModelVersionStatus(string b, bool c) {};
};
//class ModelVersionStatus: public JSON_Message {
//    public:
//        // General Vars
//        const char my_name[21] = "Model Version Status"; // Name of message for error messages
//        ModelVersionStatus(string &json_string, bool verbose_flag=false)
//            {
//                verbose = verbose_flag;
//                //JSON_Message::initialize(json_string);
//            }
//    private:
//        // Vars to read from message
//        string version;
//        string state;
//        string error_code;
//        string error_message;
//
//        // Message constants:
//        const char mmvs[21] = "model_version_status",
//                   mversion[8] = "version",
//                   mstate[6] = "state",
//                   mstatus[7] = "status",
//                   merror_code[11] = "error_code",
//                   merror_message[14] = "error_message";
//
//        // Member functions
//        void Parser() {
//            // Assert that there is a member named "model_version_status"
//            assert(document.HasMember(mmvs));
//
//            // Get the version and the state
//            version = document[mmvs][0]["version"].GetString();
//            cout << document[mmvs][0]["state"].GetString() << "\n";
//        }
//};
#endif
