#define PORT 8504
#define MODEL_NAME "full_basis"
/* This file contains utility functions for working with tensorflow-serving
 * models.
 *
 * It can query a server and check the status
 */
#include <string>
class TFServer {
  public:
      TFServer(int set_port, char *set_model_name);
      ~TFServer();
  private:
    int port;
    std::string model_name;
    //
    // This returns the URL for the model status, has information about whether
    // the model is running and if there is an error code.
    std::string get_model_status_uri();
    //
    // This returns the URI for the metadata, which has information about
    // what the name of the inputs and outputs for a given model are.
    std::string get_metadata_uri();
    /*
     * Convert a string to a character array.
     * This is not being used right now.
     */
    char* string_to_char (std::string str);
    /*
     * Initialize the interface with the tensorflow serving container
     */
    void start_interface();
    /*
     * End the interface with the tensorflow serving container
     */
    void end_interface();
    /* Send a http GET request to the curl server URI
     * Then check the resulting message to make sure the server is running
     * as expected.
     */
    void check_tfs_status();
    /* Send a http GET request to the tensorflow serving metadata URL
     */
    void get_tfs_metadata();
};

//
// C Interface
//
// We need to create c functions for each of the public c++ functions
#ifndef __cplusplus
#error You must compile with c++ compiler (g++)
#endif
extern "C" {
    class TFServer;
    // typedef aliases TFserver to TFSERVER. ..
    typedef TFServer TFSERVER;

    //Constructor
    TFSERVER* create_tfserver(int port);

    //Destructor
    void delete_tfserver(TFSERVER* tfserver);

    // The const quaificators maps from the member function to pointers in the
    // class instances
}
