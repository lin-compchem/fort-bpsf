! This pragma is the atomic number integer kind
#include "parameters.h"
program save_bp_symfuncs
   use bp_symfuncs
   use hdf5
   use h5d
   use iso_c_binding
   use h5_file_info

   implicit none
   character(len=:), allocatable :: if_path, of_path, h5_path
   real*8, allocatable :: coords(:,:,:)
   integer*ANUMKIND, allocatable :: atmnm(:,:)
   integer*NATMKIND, allocatable :: natoms(:)
   ! Get the paths from the command line
   call get_file_names(if_path, of_path)
   call read_input_file(if_path, h5_path=h5_path)
   ! Get the data from the h5 file
   call FI_init_from_ifi(h5_path) 
   call FI_get_data(coords, natoms, atmnm)
   call allocate_arrays(rad_bas, ang_bas, atmnm, mol_ids, mol2bas)
   call FI_init_outfile(of_path) 
   ! Get the cartesian gradient from the input file
   call atom_based_cartgrad(atmnm, natoms)
   
   call FI_close_input_h5()

   ! Woohoo!, the tricky parts are over now. Lets crunch some numbers
   call calculate_basis(rad_bas, ang_bas, coords, atmnm, natoms, max_atoms, &
                        num_geoms, FInum_of_els, mol_ids, mol2bas)
   call save_basis(rad_bas, ang_bas, num_els, mol_ids, mol2bas)
   call clean_up()
   return
end program save_bp_symfuncs




