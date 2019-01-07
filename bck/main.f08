module h5_file_info
   use h5d
   use hdf5
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

   type basis_pointer
   real*8, dimension(:,:), pointer :: b
   end type basis_pointer
   type(basis_pointer), dimension(FInum_els) :: ang_bas, rad_bas

   ! first dim is x/y length, second el number
   integer, allocatable :: ang_dims(:,:), rad_dims(:,:) 

   ! These are the names in the hdf5 file
   character(len=13), parameter :: atmnm_ds_name = "atomic_number"
   character(len=16), parameter :: coord_ds_name = "cartesian_coords"
   character(len=9), parameter ::  natom_ds_name = "num_atoms"

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

   subroutine allocate_basis(basis_size, pbasis)
      ! Allocate the pointers to a basis array
      implicit none
      ! I/O Vars
      real*8, dimension(:,:), pointer :: pbasis
      integer, intent(in) :: basis_size(2)

      ! Local vars
      real*8, allocatable, target :: basis(:, :)

      ! Begin
      allocate(basis(basis_size(1), basis_size(2)))
      pbasis => basis
   end subroutine allocate_basis

   subroutine create_rad_ang_bas(rad_bas, ang_bas, atmnms)
      ! Create the radial and angular basis
      use bp_symfuncs, only : num_els, els, radbas_length, angbas_length
      implicit none
      ! I/O vars
      type(basis_pointer), dimension(num_els), intent(out) :: ang_bas, rad_bas
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
         call alocate_basis(rad_bas%b, rad_dims(:,i))
      enddo
      
      ! Ang basis
      do i=1, FInum_els
         call get_num_of_els(els(i), atmnms, Finum_of_els(i))
         ang_dims(:,i) = [radbas_length(i), Finum_of_els(i)]
         call alocate_basis(ang_bas%b, ang_dims(:,i))
      enddo


   end subroutine create_rad_ang_bas
end module h5_file_info


program save_bp_symfuncs
   use bp_symfuncs
   use hdf5
   use h5d
   use iso_c_binding
   use h5_file_info

   implicit none
   character(len=*), parameter :: h5_path = './test_files/behler_3b-cartgeom.h5'
   real*8, allocatable :: coords(:,:,:)
   integer*1, allocatable :: atmnm(:,:)
   integer*2, allocatable :: natoms(:)
   integer i, j, k

   ! Get the data from the h5 file
   ! call initialize_rs_eta(0.8d0, 8.0d0, 24)
   call FI_init_from_ifi(h5_path) 
   call FI_get_data(coords, natoms, atmnm)
   call FI_clean_inputs()

   !
   call initialize_element_pars() ! must be called before allocate_basis
   call allocate_basis(rad_bas, ang_bas)
   call get_num_of_els(1, atmnm, i)
   call get_num_of_els(8, atmnm, j)
   call calculate_basis_size()
   

   return
end program save_bp_symfuncs

