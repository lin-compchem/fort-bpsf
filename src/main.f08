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
   ! call initialize_rs_eta(0.8d0, 8.0d0, 24)
   call FI_init_from_ifi(h5_path) 
   call FI_get_data(coords, natoms, atmnm)
   call FI_clean_inputs()

   ! Initialize the arrays
   !call initialize_element_pars() ! must be called before allocate_arrays
   call allocate_arrays(rad_bas, ang_bas, atmnm, mol_ids, mol2bas)

   ! Woohoo!, the tricky parts are over now. Lets crunch some numbers
   call calculate_basis(rad_bas, ang_bas, coords, atmnm, natoms, max_atoms, &
                        num_geoms, FInum_of_els, mol_ids, mol2bas)
   call save_basis(rad_bas, ang_bas, num_els, of_path, mol_ids, mol2bas)
   call deallocate_arrays()
   return
end program save_bp_symfuncs




