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
      TFServer(int set_port, char* my_model_name);
      // The below gets the energy
      void sendBPSF(double *basis, int *max_bas, int *max_atom,
              int *num_bas, int *num_atom, int *num_el, double *energy);
      ~TFServer();
  private:
    int port; // Port for Rest API of Tensorflow Server
    std::string model_name; // Name for Rest API of Tensorflow Server
    //
    // This returns the URL for the model status, has information about whether
    // the model is running and if there is an error code.
    std::string get_model_status_uri();
    //
    // This returns the URI for the metadata, which has information about
    // what the name of the inputs and outputs for a given model are.
    std::string get_metadata_uri();
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
    /* Print the basis to make sure that we are communicating the correct
     * information
     */
    void print_basis(double *basis, int *max_bas, int *max_atom,
        int *num_bas, int *num_atom, int *num_el);
    /* Print one element from the above function */
    void print_el_basis(double *basis, int max_bas, int max_atom, int num_bas,
                    int num_atom);
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
    //TFSERVER* create_tfserver(int port);
    TFSERVER* create_tfserver(int port, char* model_name);

    //Destructor
    void delete_tfserver(TFSERVER* tfserver);

    // The const quaificators maps from the member function to pointers in the
    // class instances
    void tfs_bpsf_energy(TFSERVER* tfserver,  double *basis,int *max_bas,
        int *max_atom, int *num_bas, int *num_atom, int *num_el,
        double *energy);
}
