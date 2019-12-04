#include "parameters.h"
module h5_file_info
    use h5d
    use hdf5
    use h5s
    use bp_symfuncs, only : basis, mol_id
    implicit none
    ! Contains the h5 file data
 
    ! Dimensions
    integer num_geoms ! one dimension of h5 arrays
    integer max_atoms ! another dimention of the h5 arrays
    integer, parameter :: coord_ndims = 3
    integer(HID_T) ifi, ofi ! input and output file handles
    integer(HSIZE_T) :: coord_dims(coord_ndims), coord_maxdims(coord_ndims)
    integer(HSIZE_T) :: natom_dims(1)
    integer(HSIZE_T) :: atmnm_dims(2)
    ! So an array is set up:
    ! [3, max_atoms, max_files]
    !integer :: atmn_dims(2)
    !integer :: natom_dims
    integer, parameter :: FInum_els = 2
    integer, allocatable :: FInum_of_els(:)
 
    ! The h5 output vars
    type(basis), dimension(FInum_els) :: ang_bas, rad_bas
    type(mol_id), dimension(FInum_els) :: mol_ids

    ! dimensions( 3, num_atoms, num_els)
    double precision, allocatable :: cart_grads(:,:,:)
    double precision, allocatable :: tmp_grads(:,:,:)
    ! The below is an array with dimensions(2, number_of_elements, num_geoms)
    ! The 2 are the range of indicies in the larger basis funciton array.
    ! I.e. if in the list of oxygen basis functions, molecule 6 was from indices
    ! 400 to 405 (fortran style index), then (:, 2, 6) = (399 to 405)
    ! Thats python style, would access from 0 based array terms:
    !     399 400 401 402 403 404
    integer*4, allocatable :: mol2bas(:,:,:)
 
    ! first dim is x/y length, second el number for ang/rad dims
    ! B2m dims is just vector with num of els in all geoms as the value
    ! m2b dims is the dimensions for mol2bas
    integer(kind=8), allocatable :: ang_dims(:,:), rad_dims(:,:), b2m_dims(:)
    integer(kind=8) :: radg_dims(4,FInum_els), angg_dims(4,FInum_els)
    integer(kind=8) :: m2b_dims(3)
 
    ! These are the names in the hdf5 file
    character(len=13), parameter :: atmnm_ds_name = "atomic_number"
    character(len=16), parameter :: coord_ds_name = "cartesian_coords"
    character(len=9), parameter ::  natom_ds_name = "num_atoms"

    character(len=18), dimension(2), parameter :: ofi_rad_names = ["h_radial_sym_funcs", "o_radial_sym_funcs"]
    character(len=19), dimension(2), parameter :: ofi_ang_names = ["h_angular_sym_funcs", "o_angular_sym_funcs"]
    character(len=19), dimension(2), parameter :: ofi_b2m_names = ["h_basis_to_molecule", "o_basis_to_molecule"]
    character(len=24), dimension(2), parameter :: ofi_rgrad_names = ["h_rad_cartesian_gradient", "o_rad_cartesian_gradient"]
    character(len=24), dimension(2), parameter :: ofi_agrad_names = ["h_ang_cartesian_gradient", "o_ang_cartesian_gradient"]
    character(20), dimension(2), parameter :: ofi_cgrad_names = &
    ["h_reference_gradient", "o_reference_gradient"]
    ! These are the handles for the input datasets
    integer(HID_T) atmnm_in, coord_in, natom_in, grad_in
    integer(HID_T) atomnm_out, coord_out, natom_out, grad_out
    integer(HID_T) coord_space
 
    ! These are the dimensions for the arrays
    
 
    ! Types for the dataset
    ! integer, parameter :: atmnm_type = 
    ! integer*NATMKIND num_atoms
    ! integer*ANUMKIND atm_num
    ! real*8 coords

contains
    ! Opens the input file and initializes the module level variables
subroutine FI_init_from_ifi(if_path)
    use bp_symfuncs, only: grad_key 
    implicit none
    ! IO Vars
    character(len=*), intent(in) :: if_path
    ! Local vars
    integer err
    logical file_exists
    ! Assert that the file exists
    inquire(file=if_path, exist=file_exists)
    if (.not. file_exists) goto 1050
    print *, if_path
    ! Initialize HDF5 interface
    call h5open_f(err)
    if (err .ne. 0) goto 1000
    !
    ! Open the file
    call h5fopen_f(if_path, H5F_ACC_RDONLY_F, ifi, err)
    if (err .ne. 0) goto 1010
    !
    ! open the datasets
    call h5dopen_f(ifi, atmnm_ds_name, atmnm_in, err)
    call h5dopen_f(ifi, coord_ds_name, coord_in, err)
    call h5dopen_f(ifi, natom_ds_name, natom_in, err)
    if (err .ne. 0) goto 1020
    if (allocated(grad_key)) then
        call h5dopen_f(ifi, grad_key, grad_in, err)
    endif
    if (err .ne. 0) goto 1060
    !
    ! Retrieve the dataset creation property list, and print the
    ! storage layout.
    !
    ! This is kept here for future reference but it is contiguous
    !
    ! call h5dget_create_plist_f(ds_id, dcpl, err)
    ! if (err .ne. 0) stop 'h3'
    ! call h5pget_layout_f(dcpl, layout, err)
    ! if (err .ne. 0) stop 'h4'
    ! WRITE(*,'(/,"Storage layout for ", A," is: ")', ADVANCE='NO') dataset
    ! IF(layout.EQ.H5D_COMPACT_F)THEN
    !    WRITE(*,'("H5D_COMPACT_F",/)')
    ! ELSE IF (layout.EQ.H5D_CONTIGUOUS_F)THEN
    !    WRITE(*,'("H5D_CONTIGUOUS_F",/)')
    ! ELSE IF (layout.EQ.H5D_CHUNKED_F)THEN
    !    WRITE(*,'("H5D_CHUNKED",/)')
    ! ENDIF
    
    !
    ! Get the dimensions of the dataset
    !
    call h5dget_space_f(coord_in, coord_space, err)
    if (err .ne. 0) goto 1030
    ! Can use below to find coord ndimf you didn't know
    !
    ! call h5Sget_simple_extent_ndims_f(coord_space, an_num_dims, err)
    ! allocate(an_dims(an_num_dims))
    ! allocate(max_dims(an_num_dims))
    ! max_dims = an_num_dims
    call h5Sget_simple_extent_dims_f(coord_space, coord_dims, coord_maxdims, err)
    if (err .ne. coord_ndims) goto 1040
    call h5sclose_f(coord_space, err)
    ! Parse and let the people know what the dimensions of their data
    num_geoms = coord_dims(3)
    max_atoms = coord_dims(2)
    natom_dims(1) = num_geoms
    atmnm_dims(1) = max_atoms
    atmnm_dims(2) = num_geoms
 
    print *, "From input file ", if_path
    print *, "Number of geometries: ", num_geoms
    print *, "Maximum dimension for number of atoms: ", max_atoms
    return
    ! Error handling
1000 print *, "Error opening the h5 interface"
    stop "init_from_ifi 1"
1010 print *, "Error opening the h5 file"
    stop "init_from_ifi 2"
1020 print *, "Error opening the coords, atomic numbers, and/or natom datasets in h5 file"
    stop "init_from_ifi 3"
1030 print *, "Error opening the coord space in h5 file"
    stop "init_from_ifi 4"
1040 print *, "Error getting coordinate dimensions"
    stop "init_from_ifi 5"
    stop "init_from_ifi 4"
1050 print *, "Error, geom file does not exist at path: ", if_path
    stop "init_from_ifi 6"
1060 print *, "Error opening the cartesian gradients from the file: ", if_path
    print *, "attempted input key: ", grad_key
    stop "init_from_ifi 7"
end subroutine FI_init_from_ifi
 
! This subroutine gets data from the HDF5 file and copies it to arrays in
! memory.
! 
! It also fills the FInum_of_els array and the gradient array if requested.
subroutine FI_get_data(coords, natoms, atmnms)
    implicit none
    ! IO Vars
    real*8, intent(out), allocatable :: coords(:,:,:)
    integer*NATMKIND, intent(out), allocatable :: natoms(:)
    integer*ANUMKIND, intent(out), allocatable :: atmnms(:,:)
    integer :: err
    ! Begin
    allocate(coords(3, max_atoms, num_geoms))
    allocate(atmnms(max_atoms, num_geoms))
    allocate(natoms(num_geoms))
    print *, natom_dims, atmnm_dims
    call h5dread_f(coord_in, H5T_NATIVE_DOUBLE, coords, coord_dims, err)
    if (err .ne. 0) goto 1000
    call h5dread_f(natom_in, H5T_STD_I16LE,  natoms, natom_dims, err)
    if (err .ne. 0) goto 1010
    call h5dread_f(atmnm_in, H5T_STD_I8LE,   atmnms, atmnm_dims, err)
    if (err .ne. 0) goto 1020

    call get_num_of_els(atmnms)
    return
1000 print *, "Error reading 'cartesian_coords' into memory"
     stop 'FI_get_data 1'
1010 print *, "Error reading 'num_atoms' into memory"
     stop 'FI_get_data 2'
1020 print *, "Error reading 'atom number' into memory"
     stop 'FI_get_data 3'
end subroutine FI_get_data

! Close the input H5 datasets and file
subroutine FI_close_input_h5()
    implicit none
    ! Local variables
    integer*4 err
    ! Close the input file
    call h5dclose_f(coord_in, err)
    if (err .ne.0 ) goto 5
    call h5dclose_f(natom_in, err)
    if (err .ne.0 ) goto 6
    call h5dclose_f(atmnm_in, err)
    if (err .ne.0 ) goto 7
    call h5fclose_f(ifi, err)
    if (err .ne.0 ) goto 8
    return
    return
   ! Error handling
 5 print *, "Error closing coord file"
   stop 'FI_close_input_h5 1'
 6 print *, "Error closing natom file"
   stop 'FI_close_input_h5 2'
 7 print *, "Error closing atmnm file"
   stop 'FI_close_input_h5 3'
 8 print *, "Error closing grad file"
   stop 'FI_close_input_h5 4'
end subroutine FI_close_input_h5

! Find the number of atomic symbols matching "el" in input array 
subroutine get_num_of_el(el, atmnms, num_of_els)
    implicit none
    ! IO Vars
    integer*ANUMKIND, intent(in) :: atmnms(:,:)
    integer, intent(in) :: el
    integer, intent(out) :: num_of_els 
    ! Local Vars
    integer*ANUMKIND :: el1
    integer i, j, dims(2)
    !Begin
    dims(1:2) = shape(atmnms)
    el1 = el
    num_of_els = 0
    do i=1, dims(2)
       do j=1, dims(1)
          if (atmnms(j, i) .eq. el1) num_of_els = num_of_els + 1
          if (atmnms(j, i) .eq. 0) exit
       enddo
    enddo
end subroutine get_num_of_el
 
! Fill the num_of_els data structure for all of the elements in the data file
subroutine get_num_of_els(atmnms)
    use bp_symfuncs, only: num_els, els
    integer*ANUMKIND, intent(in) :: atmnms(:,:)
    integer i 
    allocate(FInum_of_els(num_els)) 
    ! Get the number of each element so we can fill up the gradient array.
    do i=1, num_els
       call get_num_of_el(els(i), atmnms, FInum_of_els(i))
    enddo
end subroutine get_num_of_els


subroutine allocate_arrays(rad_bas, ang_bas, atmnms, mol_ids, mol2bas)
    ! Allocate the memory for the radial and angular basis functions
    ! Also allocate the memory for the mol_id variable
    use bp_symfuncs, only : num_els, radbas_length, angbas_length
    implicit none
    ! I/O vars
    type(basis), dimension(num_els), intent(out) :: ang_bas, rad_bas
    type(mol_id), dimension(num_els), intent(out) :: mol_ids
    integer*4, allocatable, intent(out) :: mol2bas(:,:,:)
    integer*ANUMKIND, intent(in) :: atmnms(:,:)
    ! Local vars
    integer i
    ! Begin Subroutine
    allocate(rad_dims(2, num_els))
    allocate(ang_dims(2, num_els))
    allocate(b2m_dims(num_els))
    ! Rad basis
    do i=1, num_els
       rad_dims(:, i) = [radbas_length(i), Finum_of_els(i)]
       radg_dims(:, i) = [3, max_atoms, radbas_length(i), Finum_of_els(i)]
       allocate(rad_bas(i)%b(rad_dims(1,i), rad_dims(2,i)))
       allocate(rad_bas(i)%g(radg_dims(1,i), radg_dims(2,i), radg_dims(3,i), radg_dims(4,i)))
       ang_dims(:,i) = [angbas_length(i), Finum_of_els(i)]
       angg_dims(:,i) = [3, max_atoms, angbas_length(i), Finum_of_els(i)]
       allocate(ang_bas(i)%b(ang_dims(1,i), ang_dims(2,i)))
       allocate(ang_bas(i)%g(angg_dims(1,i), angg_dims(2,i), angg_dims(3,i), angg_dims(4,i)))
       allocate(mol_ids(i)%bas2mol(Finum_of_els(i)))
    enddo
    b2m_dims(:) = Finum_of_els(:)
    m2b_dims(:) = [2, num_els, num_geoms]
    allocate(mol2bas(2, num_els, num_geoms))
end subroutine allocate_arrays

! Turn the cartesian gradient from the .xyz format to the basis format
! 
! We copy the values from the tmp_grad array to the array we will save
subroutine atom_based_cartgrad(atmnms, natoms)
    use bp_symfuncs, only: num_els, els, grad_key
    ! I/O variables
    integer*ANUMKIND, intent(in) :: atmnms(:,:) ! Atomic numbers
    integer*NATMKIND, intent(in) :: natoms(:)
    
    ! Local variables
    integer el, g, a ! counters for the element, geometry, and number of atoms
    integer b(num_els) ! counter that keeps track of atom-based index
    integer*4 err ! HDF5 error flag

    if (.not. allocated(grad_key)) return

    ! This array holds all of the grads as its likely too slow to make a
    ! whole bunch of individual h5 calls 
    allocate(tmp_grads(3, max_atoms, num_geoms))
    allocate(cart_grads(3, maxval(FInum_of_els), num_els))  
    ! Open the dataset
    call h5dread_f(grad_in, H5T_NATIVE_DOUBLE, tmp_grads, coord_dims, err)
    if (err .ne. 0) goto 1030
    
    b(:) = 0
    do g=1, num_geoms
        do a=1, natoms(g)
            do el=1, num_els
                if (atmnms(a, g) .eq. els(el)) then
                    b(el) = b(el) + 1
                    cart_grads(:,b(el),el) = tmp_grads(:,a,g)
                    exit
                endif 
            enddo
        enddo 
    enddo
    call h5dclose_f(grad_in, err)
    deallocate(tmp_grads)
    call save_cartgrads()
    return
    ! Error handling
1030 print *, "Error reading 'cartesian gradients' into memory"
     stop 'FI_get_data 4'
end subroutine atom_based_cartgrad
!
! Save the atom-based cartesian gradients to the output file
subroutine save_cartgrads()
    use bp_symfuncs, only: num_els
    
    ! Local variables
    integer(kind=4) :: error, el
    INTEGER(HID_T) :: dspace_id, dset_id ! Dataspace and dataset handles
    integer(kind=8) :: dims(2)
    
    do el=1, num_els
        dims(1) = 3
        dims(2) = FInum_of_els(el)
        CALL h5screate_simple_f(2, dims, dspace_id, error)
        if (error .ne. 0) goto 1000
    !
    ! Create the dataset with default properties.
    !
        CALL h5dcreate_f(ofi, ofi_cgrad_names(el), H5T_NATIVE_DOUBLE ,&
            dspace_id, dset_id, error)
        if (error .ne. 0) goto 1005
    !
    ! Write the data
    !
        CALL h5dwrite_f(dset_id, H5T_NATIVE_double,&
             cart_grads(:, :dims(2), el), dims, error)
        if (error .ne. 0) goto 1010
    !
    ! Close and release resources.
    !
        CALL h5dclose_f(dset_id , error)
        CALL h5sclose_f(dspace_id, error)
        if (error .ne. 0) goto 1015
    enddo
    return
    ! Error handling
1000 print *, "Error creating dataspace"
    stop 'save_cartgrads 1'
1005 print *, "Error creating dataset"
    stop 'save_cartgrads 2'
1010 print *, "Error writing dataset"
    stop 'save_cartgrads 3'
1015 print *, "Error closing dataset"
    stop 'save_cartgrads 4'
end subroutine save_cartgrads
!
! Closes the HDF5 interface
! Deallocates memory for arrays for this module based on program termination.
subroutine clean_up()
    implicit none
    ! Local vars
    integer i
    integer(kind=4) :: error

    ! Begin Subroutine
    do i=1, FInum_els
       deallocate(rad_bas(i)%b)
       deallocate(ang_bas(i)%b)
       deallocate(mol_ids(i)%bas2mol)
    enddo
    deallocate(mol2bas)
    if (allocated(cart_grads)) deallocate(cart_grads)
    !
    ! Close HDF5 FORTRAN interface.
    !
    CALL h5close_f(error)
    return
 
end subroutine clean_up

    ! Save the basis to the output hdf5 file
subroutine save_basis(rad_bas, ang_bas, num_els, mol_ids, mol2bas)
    use bp_symfuncs, only :  basis, mol_id
    !IO VARS
    integer, intent(in) :: num_els
    type(basis), dimension(num_els), intent(in) :: rad_bas, ang_bas
    type(mol_id), dimension(num_els), intent(in) :: mol_ids
    integer*4 :: mol2bas(2, num_els, num_geoms)
    
    ! Local variables
    integer(kind=4) :: error, i
    INTEGER(HID_T) :: dspace_id, dset_id ! Dataspace and dataset handles

    !
    ! Write the datasets.
    !
    !
    do i=1, num_els
    !
    ! RADIAL BASIS
    !
        CALL h5screate_simple_f(2, rad_dims(:,i), dspace_id, error)
        if (error .ne. 0) goto 1000
    !
    ! Create the dataset with default properties.
    !
        CALL h5dcreate_f(ofi, ofi_rad_names(i), H5T_NATIVE_DOUBLE , dspace_id, &
           dset_id, error)
        if (error .ne. 0) goto 1005
    !
    ! Write the data
    !
        CALL h5dwrite_f(dset_id, H5T_NATIVE_double, rad_bas(i)%b(:,:), rad_dims(:,i), error)
        if (error .ne. 0) goto 1010
    !
    ! Close and release resources.
    !
        CALL h5dclose_f(dset_id , error)
        CALL h5sclose_f(dspace_id, error)
        if (error .ne. 0) goto 1015
    !
    ! RADIAL Gradient
    !
       CALL h5screate_simple_f(4, radg_dims(:,i), dspace_id, error)
       if (error .ne. 0) goto 1000
    !
    ! Create the dataset with default properties.
    !
        CALL h5dcreate_f(ofi, ofi_rgrad_names(i), H5T_NATIVE_DOUBLE , dspace_id, &
           dset_id, error)
        if (error .ne. 0) goto 1005
    !
    ! Write the data
    !
        CALL h5dwrite_f(dset_id, H5T_NATIVE_double, rad_bas(i)%g(:,:,:,:), radg_dims(:,i), error)
        if (error .ne. 0) goto 1010
    !
    ! Close and release resources.
    !
        CALL h5dclose_f(dset_id , error)
        CALL h5sclose_f(dspace_id, error)
        if (error .ne. 0) goto 1015
    !
    ! ANGULAR BASIS
    !
    ! h5screate_simple_f(rank, dims, dspace_id, error)
       CALL h5screate_simple_f(2, ang_dims(:,i), dspace_id, error)
       if (error .ne. 0) goto 1020
    !
    ! Create the dataset with default properties.
    !
        CALL h5dcreate_f(ofi, ofi_ang_names(i), H5T_NATIVE_DOUBLE , dspace_id, &
           dset_id, error)
        if (error .ne. 0) goto 1025
    !
    ! Write the data
    !
        CALL h5dwrite_f(dset_id, H5T_NATIVE_double, ang_bas(i)%b(:,:), ang_dims(:,i), error)
        if (error .ne. 0) goto 1030
    !
    ! Close and release resources.
    !
        CALL h5dclose_f(dset_id , error)
        CALL h5sclose_f(dspace_id, error)
        if (error .ne. 0) goto 1035
    !
    ! ANGULAR Gradient
    !
    ! h5screate_simple_f(rank, dims, dspace_id, error)
       CALL h5screate_simple_f(4, angg_dims(:,i), dspace_id, error)
       if (error .ne. 0) goto 1000
    !
    ! Create the dataset with default properties.
    !
        CALL h5dcreate_f(ofi, ofi_agrad_names(i), H5T_NATIVE_DOUBLE , dspace_id, &
           dset_id, error)
        if (error .ne. 0) goto 1005
    !
    ! Write the data
    !
        CALL h5dwrite_f(dset_id, H5T_NATIVE_double, ang_bas(i)%g(:,:,:,:), angg_dims(:,i), error)
        if (error .ne. 0) goto 1010
    !
    ! Close and release resources.
    !
        CALL h5dclose_f(dset_id , error)
        CALL h5sclose_f(dspace_id, error)
        if (error .ne. 0) goto 1015
    !
    ! BAS2MOL
    !
    ! h5screate_simple_f(rank, dims, dspace_id, error)
       CALL h5screate_simple_f(1, b2m_dims(i), dspace_id, error)
       if (error .ne. 0) goto 1060
    !
    ! Create the dataset with default properties.
    !
        CALL h5dcreate_f(ofi, ofi_b2m_names(i), H5T_STD_I32LE , dspace_id, &
           dset_id, error)
        if (error .ne. 0) goto 1065
    !
    ! Write the data
    !
        CALL h5dwrite_f(dset_id, H5T_STD_I32LE, [mol_ids(i)%bas2mol(:)], &
            [b2m_dims(i)], error)
        if (error .ne. 0) goto 1070
    !
    ! Close and release resources.
    !
        CALL h5dclose_f(dset_id , error)
        CALL h5sclose_f(dspace_id, error)
        if (error .ne. 0) goto 1075
    enddo
    !
    ! Write the mol2bas data structure
    !
    CALL h5screate_simple_f(3, m2b_dims, dspace_id, error)
    if (error .ne. 0) goto 1040
    !
    ! Create the dataset with default properties.
    !
    CALL h5dcreate_f(ofi, 'molecule_to_basis_inds', H5T_STD_I32LE , &
       dspace_id, dset_id,  error)
    if (error .ne. 0) goto 1045
    !
    ! Write the data
    !
    CALL h5dwrite_f(dset_id, H5T_STD_I32LE, mol2bas(:,:,:), m2b_dims, error)
    if (error .ne. 0) goto 1050
    !
    ! Close and release resources.
    !
    CALL h5dclose_f(dset_id , error)
    CALL h5sclose_f(dspace_id, error)
    if (error .ne. 0) goto 1055
    !
    ! Close the file.
    !
    CALL h5fclose_f(ofi, error)
    return
1000 print *, "Error writing dataset"
    stop 'save_basis 1'
1005 print *, "Error writing dataset"
    stop 'save_basis 2'
1010 print *, "Error writing dataset"
    stop 'save_basis 3'
1015 print *, "Error writing dataset"
    stop 'save_basis 4'
1020 print *, "Error writing dataset"
    stop 'save_basis 5'
1025 print *, "Error writing dataset"
    stop 'save_basis 6'
1030 print *, "Error writing dataset"
    stop 'save_basis 7'
1035 print *, "Error writing dataset"
    stop 'save_basis 8'
1040 print *, "Error writing mol2bas dataset"
    stop 'save_basis 9'
1045 print *, "Error writing mol2bas dataset"
    stop 'save_basis 10'
1050 print *, "Error writing mol2bas dataset"
    stop 'save_basis 11'
1055 print *, "Error writing mol2bas dataset"
    stop 'save_basis 12'
1060 print *, "Error writing b2m dataset"
    stop 'save_basis 13'
1065 print *, "Error writing b2m dataset"
    stop 'save_basis 14'
1070 print *, "Error writing b2m dataset"
    stop 'save_basis 15'
1075 print *, "Error writing b2m dataset"
    stop 'save_basis 16'
end subroutine save_basis

! Gets the file names from the command line. The input file path should be
! the first argument and the output file path should be the second.
subroutine get_file_names(h5_path, of_path)
   implicit none
   !I/O Vars
   character(len=:),allocatable, intent(inout) :: h5_path, of_path
   character(len=400) :: arg
   integer :: length, status

   !Begin
   call get_command_argument(1, value=arg, length=length, status=status)
   if (status .gt. 0) then
      print *, "Please specify the input file as the first argument"
      stop 'get_file_names 1'
   elseif (status .lt. 0) then
      print *, "Error retrieving input file command line argument"
      stop 'get_file_names 2'
   else
      h5_path = arg(:length)
   endif
   call get_command_argument(2, value=arg, length=length, status=status)
   if (status .gt. 0) then
      print *, "Please specify the output file as the second argument"
      stop 'get_file_names 3'
   elseif (status .lt. 0) then
      print *, "Error retrieving output command line argument"
      stop 'get_file_names 4'
   else
      of_path = arg(:length)
   endif

   return
end subroutine get_file_names
! Opens the output HDF5 file
subroutine FI_init_outfile(of_path)
    character(len=*), intent(in) :: of_path

    ! Local variables
    integer(kind=4) :: error
    
    !
    ! Open a new file.
    !
    call h5fcreate_f(of_path, H5F_ACC_TRUNC_F, ofi, error)
    if (error .ne. 0)  goto 2000
    return

    ! Error handling
    2000 print *, "error opening output file: ", of_path
    stop 'save_basis 17'
end subroutine FI_init_outfile
end module h5_file_info
