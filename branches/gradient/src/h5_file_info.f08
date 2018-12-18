module h5_file_info
    use h5d
    use hdf5
    use h5s
    use bp_symfuncs, only : basis
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
 
    type(basis), dimension(FInum_els) :: ang_bas, rad_bas
 
    ! first dim is x/y length, second el number
    integer(kind=8), allocatable :: ang_dims(:,:), rad_dims(:,:)
 
    ! These are the names in the hdf5 file
    character(len=13), parameter :: atmnm_ds_name = "atomic_number"
    character(len=16), parameter :: coord_ds_name = "cartesian_coords"
    character(len=9), parameter ::  natom_ds_name = "num_atoms"

    character(len=18), dimension(2), parameter :: ofi_rad_names = ["h_radial_sym_funcs", "o_radial_sym_funcs"]
    character(len=19), dimension(2), parameter :: ofi_ang_names = ["h_angular_sym_funcs", "o_angular_sym_funcs"]

 
    ! These are the handles for the input datasets
    integer(HID_T) atmnm_in, coord_in, natom_in
    integer(HID_T) atomnm_out, coord_out, natom_out
    integer(HID_T) coord_space
 
    ! These are the dimensions for the arrays
    
 
    ! Types for the dataset
    ! integer, parameter :: atmnm_type = 
    ! integer*2 num_atoms
    ! integer*1 atm_num
    ! real*8 coords

 contains
    subroutine FI_init_from_ifi(if_path)
       ! This subroutine opens the input file and initializes the module level variables
       implicit none
       ! IO Vars
       character(len=*), intent(in) :: if_path
       ! Local vars
       integer err
       ! Begin subroutine!
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
       !call h5dget_space_f(atmnm_in, atmnm_space, err)
       !call h5dget_space_f(natom_in, natom_space, err)
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
 1000  print *, "Error opening the h5 interface"
       stop "init_from_ifi 1"
 1010  print *, "Error opening the h5 file"
       stop "init_from_ifi 2"
 1020  print *, "Error opening the coords, atomic numbers, and/or natom datasets in h5 file"
       stop "init_from_ifi 3"
 1030  print *, "Error opening the coord space in h5 file"
       stop "init_from_ifi 4"
 1040  print *, "Error getting coordinate dimensions"
       stop "init_from_ifi 5"
 
    end subroutine FI_init_from_ifi
 
    subroutine FI_get_data(coords, natoms, atmnms)
       ! This subroutine just copies the data from the h5py file to
       ! arrays
       implicit none
       ! IO Vars
       real*8, intent(out), allocatable :: coords(:,:,:)
       integer*2, intent(out), allocatable :: natoms(:)
       integer*1, intent(out), allocatable :: atmnms(:,:)
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
       return
  1000 print *, "Error reading 'cartesian_coords' into memory"
       stop 'FI_get_data 1'
  1010 print *, "Error reading 'num_atoms' into memory"
       stop 'FI_get_data 2'
  1020 print *, "Error reading 'atom number'b into memory"
       stop 'FI_get_data 3'
    end subroutine FI_get_data
 
    subroutine FI_clean_inputs()
       implicit none
    ! Clean up loose ends with the h5py file
       ! Local variables
       integer*4 err
       ! Close the input file
       call h5dclose_f(coord_in, err)
       call h5dclose_f(natom_in, err)
       call h5dclose_f(atmnm_in, err)
       call h5fclose_f(ifi, err)
    end subroutine FI_clean_inputs
 
    ! function FI_find_num_atoms_of_el(el)
    !    ! Find the number of occurences of the atomic number in the array
    !    implicit none
    !    ! IO Vars
    !    integer*1 el
    !    integer FI_find_num_atoms_of_el
 
    !    ! Local vars
    !    integer i, j
    !    integer*1 rbuf
 
    !    !Begin
    !    do i=1, num_geoms
    !       do j=1, max_atoms
    !           call h5dread_f(atmnm_in, H5T_NATIVE_SHORT, atmnm_space)
    !       enddo
    !    enddo
 
 
    ! end function FI_find_num_atoms_of_el
    subroutine get_num_of_els(el, atmnms, num_of_els)
       ! Find the number of atomic symbols matching "el" in input array
       implicit none
       ! IO Vars
       integer*1, intent(in) :: atmnms(:,:)
       integer, intent(in) :: el
       integer, intent(out) :: num_of_els 
       ! Local Vars
       integer*1 :: el1
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
    end subroutine get_num_of_els
 
 
    subroutine create_rad_ang_basis(rad_bas, ang_bas, atmnms)
       ! Create the radial and angular basis
       use bp_symfuncs, only : num_els, els, radbas_length, angbas_length
       implicit none
       ! I/O vars
       type(basis), dimension(num_els), intent(out) :: ang_bas, rad_bas
       integer*1, intent(in) :: atmnms(:,:)
       ! Local vars
       integer i
       ! Begin Subroutine
       ! FInum_els = num_els
       allocate(FInum_of_els(FInum_els))
       allocate(rad_dims(2, FInum_els))
       allocate(ang_dims(2, FInum_els))
       ! Rad basis
       do i=1, FInum_els
          call get_num_of_els(els(i), atmnms, Finum_of_els(i))
          rad_dims(:, i) = [radbas_length(i), Finum_of_els(i)]
          allocate(rad_bas(i)%b(rad_dims(1,i), rad_dims(2,i)))
          ang_dims(:,i) = [angbas_length(i), Finum_of_els(i)]
          allocate(ang_bas(i)%b(ang_dims(1,i), ang_dims(2,i)))
       enddo
 
 
    end subroutine create_rad_ang_basis


    subroutine save_basis(rad_bas, ang_bas, num_els, of_path)
    ! Save the basis to the output hdf5 file
      use bp_symfuncs, only :  basis
      implicit none
      character(len=*), intent(in) :: of_path
      integer :: num_els
      type(basis), dimension(num_els) :: rad_bas, ang_bas
      
      ! Local variables
      integer(kind=4) :: error, i
      INTEGER(HID_T) :: dspace_id, dset_id ! Dataspace and dataset handles
   !
   ! Initialize FORTRAN interface.
   !
      call h5open_f(error)
   !
   ! Open a new file.
   !
      call h5fcreate_f(of_path, H5F_ACC_TRUNC_F, ofi, error)
   !
   ! Write the datasets.
   !
   !
      do i=1, num_els
      !
      ! RADIAL BASIS
      !
      ! h5screate_simple_f(rank, dims, dspace_id, error)
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
          if (error .ne. 0) goto 1015
      enddo
  !
  ! Close the file.
  !
  CALL h5fclose_f(ofi, error)
  !
  ! Close FORTRAN interface.
  !
   CALL h5close_f(error)
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

    end subroutine save_basis

subroutine get_file_names(h5_path, of_path)
   ! Gets the file names from the command line. The input file path should be
   ! the first argument and the output file path should be the second.
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
 end module h5_file_info
