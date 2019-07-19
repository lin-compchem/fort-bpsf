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

#ifndef TFS_JSON
#define TFS_JSON
class JSON_Message {
    public:
        // General Vars
        rapidjson::Document document; // This is the object to parse the JSON string into
        const char my_name[8] = "generic"; // Name of message for error messages
        int parse_status = -1;
        // Member functions
        JSON_Message();
        JSON_Message(std::string &json_string, bool verbose_flag=false);
    protected:
        bool verbose = false; // Flag for verbosity
        void setupDocument(std::string &json_string); 
        virtual int parseDocument();
        void parseFail();
        void initialize(std::string &json_string);
};
#endif
#ifndef TFS_MVS
#define TFS_MVS
class ModelVersionStatus: public JSON_Message {
    public:
        // General Vars
        const char my_name[21] = "Model Version Status"; // Name of message for error messages
        ModelVersionStatus(std::string &json_string, bool verbose_flag=false);
        void checkStatus();
    private:
        // Vars to read from message
        std::string version;
        std::string state;
        std::string error_code;
        std::string error_message;

        // Message constants:
        const char mmvs[21] = "model_version_status",
                   mversion[8] = "version",
                   mstate[6] = "state",
                   mstatus[7] = "status",
                   merror_code[11] = "error_code",
                   merror_message[14] = "error_message",
                   good_code[3] = "OK",
                   good_state[10] = "AVAILABLE";

        // Member functions
        int parseDocument(); 
        void printStatus();
};
#endif
