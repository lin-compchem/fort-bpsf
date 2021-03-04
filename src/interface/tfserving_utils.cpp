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
#include <fstream>
#include <stdlib.h>
// Program includes
#include "tfs_messages.hpp"
#include "tfserving_utils.hpp"
#include "curler.hpp"
// rapidjson
#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/prettywriter.h>
using namespace std;
using namespace rapidjson;
/*
 * Constructor
 */
TFServer::TFServer(int set_port, char* my_model_name):port(set_port) {
   model_name = my_model_name;
   start_interface();  
}
/*
 * Destructor
 */
TFServer::~TFServer() {
    end_interface();
}
/*
 * This returns the URL for the model status, has information about whether
 * the model is running and if there is an error code.
 */
string TFServer::get_model_status_uri() {
    string uri = "http://localhost:";
    uri += to_string(port);
    uri += "/v1/models/" + model_name;
    return uri;
}
/*
 * This returns the URL for the model status, has information about whether
 * the model is running and if there is an error code.
 */
string TFServer::get_predict_uri() {
    string uri = "http://localhost:";
    uri += to_string(port);
    uri += "/v1/models/" + model_name;
    uri += ":predict";
    return uri;
}
/*
 * This returns the URI for the metadata, which has information about
 * what the name of the inputs and outputs for a given model are.
 */
string TFServer::get_metadata_uri() {
    string uri = "http://localhost:";
    uri += to_string(port);
    uri += "/v1/models/" + model_name + "/metadata";
    return uri;
}
/*
 * Send a http GET request to the curl server URI
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
/* 
 * Send a http GET request to the tensorflow serving metadata URL
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
}
/*
 * This accepts the basis arrays and gets the energy and gradients
 * from the tensorflow server
 */
void TFServer::sendBPSF(double *basis, int *max_bas, int *max_atom,
                   int *num_bas, int *num_atom, int *num_el, double *energy) {
    if (verbose) { cout << "Entering sendBPSF" << endl;}
    // Write the JSON string to outstr 
    string basis_json = serialize_bpsf(basis, max_bas, max_atom, num_bas,
            num_atom, num_el);

    // Send the string to the server and get the answer
    string post_str = send_basis(basis_json.c_str());

    // Parse the response with the messenger object
    EnGradMessage engrad(post_str, true);
    energy[0] = engrad.parseEnergy();
    
    if (verbose) {cout << endl << "Leaving sendBPSF" << endl;}
    return;
}
/*
 * This accepts the basis arrays and gets the energy and gradients
 * from the tensorflow server
 */
void TFServer::sendBPSF(double *basis, int *max_bas, int *max_atom,
                   int *num_bas, int *num_atom, int *num_el, double *energy,
                   double *gradient) {
    // if (verbose) { cout << "Entering sendBPSF" << endl;}
    // Write the JSON string to outstr 
    string basis_json = serialize_bpsf(basis, max_bas, max_atom, num_bas,
            num_atom, num_el);

    // Send the string to the server and get the answer
    string post_str = send_basis(basis_json.c_str());

    // Parse the response with the messenger object
    EnGradMessage engrad(post_str, true);
    energy[0] = engrad.parseEnergy();
    // Gradient
    string gnames[2] = {"h_basis_grad", "o_basis_grad"};
    for (int i=0; i < num_el[0]; ++i) {
        int s = i * *max_bas * max_atom[0];
        for(int j=0; j < num_atom[i]; ++j) { 
            int b = j * max_bas[0] + s;
            for (int k=0; k < num_bas[i]; ++k) {
                int idx = b + k;
               gradient[idx] = engrad.document["outputs"][gnames[i].c_str()][j][k].GetDouble(); 
            }
        }
    }
    
    if (verbose) {cout << endl << "Leaving sendBPSF" << endl;}
    return;
}
/*
 * Take the basis arrays and serialize them into  JSON STring
 */
string TFServer::serialize_bpsf(double *basis, int *max_bas, int *max_atom,
                   int *num_bas, int *num_atom, int *num_el) {
    StringBuffer outstr; // This holds the JSON string while the writer writes
    Writer<StringBuffer> writer(outstr); // Wites JSON String
   
    // Inititialize basis keys
    string basis_keys[2] = {"h_basis", "o_basis"};

    // This is the main object
    writer.StartObject();
    writer.Key("inputs"); // Object for list of input tensors
    writer.StartObject(); // Everything below is corresponds to 'inputs' key
 
    // This is the basis object
    for (int i=0; i < num_el[0]; ++i) {
        int j = i * *max_bas * *max_atom;
        writer.Key(basis_keys[i].c_str());
        writer.StartArray(); // Begin of array for key basis
        write_el_basis(&basis[j], max_bas[0], max_atom[0], num_bas[i],
               num_atom[i], writer);
        writer.EndArray(); // end of array for key basis
    }
    // Write the bas2mol keys
    write_bas2mol(num_el[0], num_atom, writer);

    //Sanity check for debugging

    writer.EndObject(); // End inputs object
    writer.EndObject(); // End main object
    
    // For debugging
    /*
    if (TFServer::verbose) {
        cout << "Serialized basis:" << outstr.GetString() << endl;
    }
    */

    return outstr.GetString();
}
//
// This tesets the input structure for the bpnn with fake data.
void TFServer::ModelTest1() {
    StringBuffer outstr; // This holds the JSON string while the writer writes
    Writer<StringBuffer> writer(outstr); // Wites JSON String
    string basis_keys[2] = {"h_basis", "o_basis"};
    string m2b_keys[2] = {"h_bas2mol", "o_bas2mol"}; 
    string url = get_predict_uri(); 
    cout << endl << endl << "Performing Model Test 1" << endl;
    
    // This is the main object 
    writer.StartObject();
    writer.Key("inputs");
    //writer.Key("input");
    cout << "Created Main Object" << endl;
    writer.StartObject();
    // Write the basis to the object
    for (int i=0; i<2; ++i) {
        writer.Key(basis_keys[i].c_str());
        writer.StartArray();
        for (int j=0; j<2; ++j) {
            writer.StartArray();
            for (int k=1; k<4; ++k) {
                writer.Double(3*j + k);
            } 
            writer.EndArray();
        } 
        writer.EndArray();
        cout << outstr.GetString() << endl;
    } 
    // Write the mol2bas keys
    for (int i=0; i<2; ++i) {
        writer.Key(m2b_keys[i].c_str());
        writer.StartArray();
        writer.Int(1);
        writer.Int(2);
        writer.EndArray();
    }
    // Write the correction energies
    writer.EndObject();
    writer.EndObject();
   
    // End of main object 
    //writer.EndObject();
    cout << "Final Prediction String:" << endl;
    cout << outstr.GetString() << endl;
    cout << endl << endl;
    cout << "Predict URL:" << endl << url << endl;
    cout << "Performing Prediction:";

    Curler post_curl;
    post_curl.setURL(url.c_str());
    post_curl.httpPost(outstr.GetString());
   
    
    cout << "Results of POST:" << endl;
    cout << post_curl.curl_buffer;

    cout << "Time to validate results:" << endl;
    Document document;
    document.Parse(post_curl.curl_buffer.c_str());

    if (!document.HasMember("outputs")) {
        cerr << "Error, cannot read 'outputs' field in POST response" << endl;
        exit(EXIT_FAILURE);
    }
    double vals[3];
    for(int i=0; i<3; ++i) {
        vals[i] = document["outputs"][i][0].GetDouble();
    }
    cout << "Parsed Results:" << endl;
    printf("%f   %f    %f\n", vals[0], vals[1], vals[2]);
    cout << "Checking Results:" << endl;
    
    double answers[3] = {0.0, 12.0, 30.0};
    for(int i=0; i<3; ++i) {
        printf("assert val %f == %f\n", vals[i], answers[i]);
        assert(vals[i] == answers[i]);
    }

    cout << "TEST COMPLETED!" << endl;
    cout << "Ending Model Test 1" << endl;

    exit(EXIT_FAILURE);
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
// Prints the basis
// There are a lot of parameters here
// The basis object has [Fortran style] indices of
//  (max_bas, max_atom, total_els). However, many of these will be blank.
// Therefore we need to only index over the number of things that we have.
// This is done by the counters num_bas, numatom, and num_el
//
// double *basis: array of size (max_bas, max_atom, num_el)
//                this is the basis to print
void TFServer::print_basis(double *basis, int *max_bas, int *max_atom,
        int *num_bas, int *num_atom, int *num_el){
    int i, j;
    for (i=0; i < *num_el; ++i)
    {
        j = i * *max_bas * *max_atom;
        print_el_basis(&basis[j], max_bas[0], max_atom[0], num_bas[i],
               num_atom[i]);
    }
}
void TFServer::print_el_basis(double *basis, int max_bas, int max_atom,
        int num_bas, int num_atom) {
    int i, j, idx;
    int basis_counter = 0;
    for (i=0; i < num_atom; ++i) {
        basis_counter = i * max_bas;
        for (j=0; j < num_bas; ++j){
            idx = basis_counter + j;
            printf("%f           %d \n", basis[idx], idx);
        }
    }
}
//
// Write the basis set to the writer for a given element
void TFServer::write_el_basis(double *basis, int max_bas, int max_atom,
        int num_bas, int num_atom, Writer<StringBuffer> &writer) {
    int i, j, idx;
    int basis_counter = 0;
    for (i=0; i < num_atom; ++i) {
        basis_counter = i * max_bas;
        writer.StartArray();
        for (j=0; j < num_bas; ++j){
            idx = basis_counter + j;
            writer.Double(basis[idx]);
        }
        writer.EndArray();
    }
}
//
// Input:
// 
void TFServer::write_bas2mol(int num_el, int *num_atom, Writer<StringBuffer> &writer) {
//void TFServer::write_bas2mol(int num_el, int *num_atom, Writer<GenericStringBuffer<UTF8<char>, CrtAllocator>, UTF8<char>, UTF8<char>, rapidjson::CrtAllocator, 0u>& writer) {
    string b2m_keys[2] = {"h_bas2mol", "o_bas2mol"}; 
    for (int el=0; el < num_el; ++el) {
        writer.Key(b2m_keys[el].c_str());
        writer.StartArray();
        for (int atom=0; atom < num_atom[el]; ++atom) {
           writer.Int(0); 
        }
        writer.EndArray();     
    } 
    return;
}
//
// Send and arbitrary outstring to the TensorflowServing server 
string TFServer::send_basis(const char *json) {
    string url = get_predict_uri(); 
    if (verbose) {
        cout << "Performing Prediction" << endl;
        cout << "Prediction URL: " << url << endl;
    }

    /// Send the json string to the server
    Curler post_curl;
    post_curl.setURL(url.c_str());
    post_curl.httpPost(json);
   
    // Print out the response 
    if (verbose) {
        cout << "Results of POST:" << endl;
        cout << post_curl.curl_buffer;
    }

    return post_curl.curl_buffer;
}

