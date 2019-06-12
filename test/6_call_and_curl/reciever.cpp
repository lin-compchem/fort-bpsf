// Defines size_t, FILE, and fpos_t
#include <stdio.h>
// This is where size_t is defined
#include <stdlib.h>
#include <curl/curl.h>
#include <iostream>
#include <string>
#include <rapidjson/rapidjson.h>
#include <rapidjson/document.h>
#include <rapidjson/prettywriter.h>

using namespace rapidjson;
using namespace std;
// Prints a 2d array
// Note that one of the passed integers is int*2 while one is int*4
//
#ifdef __cplusplus
    extern"C" {
#endif        
void print2d_(double *arr, short *nx, int *ny);

void json_coords_(double *arr, short *nx, int *ny);

void print_basis_(double *basis, int *max_bas, int *max_atom,
        int *num_bas, int *num_atom, int *num_el);
#ifdef __cplusplus
    }
#endif 
void print_el_basis_(double *basis, int max_bas, int max_atom, int num_bas,
        int num_atom); 
//
//
//
//
//
void print_el_basis_(double *basis, int max_bas, int max_atom, int num_bas,
        int num_atom) {
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

void 
print2d_(double *arr, short *nx, int *ny) {
    int i, j;
    for(j = 0; j < *ny; ++j)
        for(i = 0; i < *nx; ++i)
        {
            int idx = j * *nx + i;
            printf("arr[%d] = %f\n", idx, arr[idx]);
        }
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
void print_basis_(double *basis, int *max_bas, int *max_atom,
        int *num_bas, int *num_atom, int *num_el){
    int i, j;
    for (i=0; i < *num_el; ++i)
    {
        j = i * *max_bas * *max_atom;
        print_el_basis_(&basis[j], max_bas[0], max_atom[0], num_bas[i],
               num_atom[i]);
    }
}

void json_coords_(double *arr, short *nx, int *ny) {
    int i, j;
    //Create a json document
    Document json_tree;
    //string name = "coordinates";
    char name[] = "coordinates";
    StringBuffer outstr;
    Writer<StringBuffer> writer(outstr);

    writer.StartObject();
    writer.Key(name);
    writer.StartArray();

    for(j = 0; j < *ny; ++j)
        for(i = 0; i < *nx; ++i)
        {
            int idx = j * *nx + i;
            printf("%f aaaaaaa\n", arr[idx]);
            writer.Double(arr[idx]);
        }
    writer.EndArray();
    writer.EndObject();

    cout << outstr.GetString() << endl;
}

/*
void json_save2d_(Document &document, string name, double *arr, int x, int y)
{
   int 
}
*/
