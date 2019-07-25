#pragma once

#include <string>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
// RapidJSON Headers
#include <rapidjson/document.h> // includes rapidjson/rapidjson.h 
#include <rapidjson/rapidjson.h>

#define PRETTYWRITER
#ifdef PRETTYWRITER
  #include <rapidjson/prettywriter.h>
#else
  #include <rapidjson/writer.h>
#endif

#ifndef TFS_JSON
#define TFS_JSON
class JSON_Message {
    public:
        // General Vars
        rapidjson::Document document; // This is the object to parse the JSON string into
        const char my_name[8] = "generic"; // Name of message for error messages
        int parse_status = -1; // Error flag
        // Member functions
        JSON_Message();
        JSON_Message(std::string &json_string, bool verbose_flag=false);
    protected:
        bool verbose = false; // Flag for verbosity
        /*
         * READ JSON VARIABLES
         *  The following subroutines are for putting JSON variables encoded
         * in &json_string into a document and then parsing the variables
         * in said document
         */
        // Have the document object parse the &json_string
        void setupDocument(std::string &json_string); 
        // Get the data out of the document object after setupDocument is called
        virtual int parseDocument();
        // Error handling for above
        void parseFail();
        // setup and parse the documents
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
//
//
class EnGradMessage: public JSON_Message {
    public:
        // General Vars
        const char my_name[24] = "Energy/Gradient Results"; // Name of message for error messages
        EnGradMessage(std::string &json_string, bool verbose_flag=false);
        double parseEnergy();
        void checkOutput(const char *name);
    private:
        std::string message; // The message that this class is initialized with

        // Message constants:
        const char mout[8] = "outputs";
        std::string mener = "correction_energies",
                    mhgrad = "h_basis_grad",
                    mograd = "o_basis_grad";

        // Member functions
        int parseDocument(); 
        void printMessage();
        void parseDoubles();
};
#endif
