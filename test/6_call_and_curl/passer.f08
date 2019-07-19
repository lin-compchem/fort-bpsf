!> This a program that is used to call the symmetry function subroutines as
!! an external library. It uses the coordinates from one of the other test
!! runs
#define MAXATOM 151
#define MAXBAS 150
#define NUMELS 2
#define NUMH 6
#define NUMO 3
#define RAD 48
#define HANG 36
#define OANG 54
program caller
    use bp_symfuncs
    implicit none
    !
    ! Setup Variables
    !
    integer*2, parameter :: natm = 9
    real*8 :: coords(3,natm)
    integer*1 :: atmnum(natm)
    character(len=:), allocatable :: in_path

    in_path = 'bp.inp'

    ! Read the .xyz coord file
    call read_xyzfile(natm, coords, atmnum)
    call print_coords(natm, coords, atmnum)

    ! This needs to always be called first to initialize the module
    call read_input_file(in_path)
    ! Now we can define variables based on how it was initialized
    call stage_calculation(natm, coords, atmnum)

end program caller
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine send_to_c(natm, coords)
    implicit none
    real*8 :: coords(3,natm)
    integer*2 :: natm
    external print_coords

    call print2d(coords, natm, 3)


end subroutine send_to_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine stage_calculation(natm, coords, atmnum)
    ! After initialize_element_pars has been called, we can
    ! go and do the calculation
    use bp_symfuncs, only : calc_bp
    use iso_fortran_env
    implicit none
    !
    ! IO Variables
    !
    integer*2 :: natm
    real*8 :: coords(3,natm)
    integer*1 :: atmnum(natm)
    !***********************************************
    ! BP Variables
    !***********************************************
    ! The number of geometries for the calculation
    integer, parameter :: num_geoms = 1
    ! The maximum number of atoms in the calculation
    integer :: max_atoms
    ! The number of each element. Think they are defined
    ! in order as H,O. DIM=(num_els)
    integer, parameter :: num_of_els(NUMELS) = [6,3]
    ! This is vestigial and can be removed to make it easeir to work
    ! with but returns the number of elements
    integer :: g_num_of_els(2) = [0,0]
    real*8, dimension(MAXBAS, MAXATOM, NUMELS) :: rad_bas, ang_bas
    real*8, dimension(3, MAXATOM, MAXBAS, MAXATOM, NUMELS) :: rad_grad, ang_grad
    !***********************************************
    ! Local variables
    !
    integer :: num_rads(2) = [RAD, RAD], &
               num_angs(2) = [HANG, OANG]
    max_atoms = MAXATOM
    print *, 'hello'
    flush(output_unit)
    ! Calculate the basis function    
    call calc_bp (natm, coords, atmnum, rad_bas, ang_bas, &
             rad_grad, ang_grad, max_atoms, g_num_of_els)

    ! Print the coords
    call send_to_c(natm, coords)
    call print_basis(rad_bas, MAXBAS, MAXATOM, num_rads, g_num_of_els, NUMELS)
    call json_coords(coords, natm, 3)
    return

end subroutine stage_calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_basis(rad_bas, ang_bas, rad_grad, ang_grad)
    use hdf5
    use h5s
    use h5d
    implicit none

    character(len=*), parameter :: out_path = 'output.h5'

    integer(hid_t) :: file_id
    integer(hid_t) :: dset_id
    integer(hid_t) :: dspace_id

    integer(kind=4) :: error

    real*8, dimension(MAXBAS, MAXATOM, NUMELS) :: rad_bas, ang_bas
    real*8, dimension(3, MAXATOM, MAXBAS, MAXATOM, NUMELS) :: rad_grad, ang_grad

    integer(kind=8), dimension(2) :: hrb_dims = [RAD, NUMH]
    integer(kind=8), dimension(2) :: orb_dims = [RAD, NUMO]

    integer(kind=8), dimension(2) :: hab_dims = [HANG, NUMH]
    integer(kind=8), dimension(2) :: oab_dims = [OANG, NUMO]

    integer(kind=8), dimension(4) :: hrg_dims = [3, MAXATOM, RAD, NUMH]
    integer(kind=8), dimension(4) :: org_dims = [3, MAXATOM, RAD, NUMO]

    integer(kind=8), dimension(4) :: hag_dims = [3, MAXATOM, HANG, NUMH]
    integer(kind=8), dimension(4) :: oag_dims = [3, MAXATOM, OANG, NUMO]
  !
  ! Initialize FORTRAN interface.
  !
  CALL h5open_f(error)

  !
  ! Create a new file using default properties.
  !
  CALL h5fcreate_f(out_path, H5F_ACC_TRUNC_F, file_id, error)

  ! *************** H RAD BASIS ********************
  CALL h5screate_simple_f(2, hrb_dims, dspace_id, error)
  CALL h5dcreate_f(file_id, "h_radial_sym_funcs", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_double, rad_bas(1:RAD,:NUMH,1), hrb_dims, error)
  CALL h5dclose_f(dset_id , error)
  CALL h5sclose_f(dspace_id, error)
  ! *************** O RAD BASIS ********************
  CALL h5screate_simple_f(2, orb_dims, dspace_id, error)
  CALL h5dcreate_f(file_id, "o_radial_sym_funcs", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_double, rad_bas(:RAD,:NUMO,2), orb_dims, error)
  CALL h5dclose_f(dset_id , error)
  CALL h5sclose_f(dspace_id, error)
  ! *************** H ANG BASIS ********************
  CALL h5screate_simple_f(2, hab_dims, dspace_id, error)
  CALL h5dcreate_f(file_id, "h_angular_sym_funcs", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_double, ang_bas(:HANG,:NUMH,1), hab_dims, error)
  CALL h5dclose_f(dset_id , error)
  CALL h5sclose_f(dspace_id, error)
  ! *************** O ANG BASIS ********************
  CALL h5screate_simple_f(2, oab_dims, dspace_id, error)
  CALL h5dcreate_f(file_id, "o_angular_sym_funcs", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_double, ang_bas(:OANG,:NUMO,2), oab_dims, error)
  CALL h5dclose_f(dset_id , error)
  CALL h5sclose_f(dspace_id, error)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                GRADIENTS
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! *************** H RAD GRAD *********************
  CALL h5screate_simple_f(4, hrg_dims, dspace_id, error)
  CALL h5dcreate_f(file_id, "h_rad_cartesian_gradient", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_double, rad_grad(:,:MAXATOM,:RAD,:NUMH,1), hrg_dims, error)
  if (error .ne. 0) goto 5000
  CALL h5dclose_f(dset_id , error)
  if (error .ne. 0) goto 5000
  CALL h5sclose_f(dspace_id, error)
  if (error .ne. 0) goto 5000
  ! *************** O RAD GRAD *********************
  CALL h5screate_simple_f(4, org_dims, dspace_id, error)
  CALL h5dcreate_f(file_id, "o_rad_cartesian_gradient", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_double, rad_grad(:,:MAXATOM,:RAD,:NUMO,2), org_dims, error)
  CALL h5dclose_f(dset_id , error)
  CALL h5sclose_f(dspace_id, error)
  ! *************** H ANG GRAD *********************
  CALL h5screate_simple_f(4, hag_dims, dspace_id, error)
  CALL h5dcreate_f(file_id, "h_ang_cartesian_gradient", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_double, ang_grad(:,:MAXATOM,:HANG,:NUMH,1), hag_dims, error)
  CALL h5dclose_f(dset_id , error)
  CALL h5sclose_f(dspace_id, error)
  ! *************** O ANG GRAD *********************
  CALL h5screate_simple_f(4, oag_dims, dspace_id, error)
  CALL h5dcreate_f(file_id, "o_ang_cartesian_gradient", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_double, ang_grad(:,:MAXATOM,:OANG,:NUMO,2), oag_dims, error)
  CALL h5dclose_f(dset_id , error)
  CALL h5sclose_f(dspace_id, error)

  !
  ! Close the file.
  !
  CALL h5fclose_f(file_id, error)

  !
  ! Close FORTRAN interface.
  !
  CALL h5close_f(error)

  return
5000 print *, 'uhoh'
stop 'i shouldnt be here'

end subroutine write_basis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine print_coords(natm, coords, atmnum)
    ! Print the coordinates to stdout in xyz file format
    implicit none
    integer*2 :: natm
    real*8 :: coords(3,natm)
    integer*1 :: atmnum(natm)

    integer :: i

    print *, "Input coordinates:"
    print *, natm
    print *, '\n'
    do i=1, natm
       print 10, atmnum(i), coords(:,i)
 10 format(I2, 5x, 3(F12.8,5x))
    enddo
end subroutine print_coords
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_xyzfile(natm, coords, atmnum)
    ! Read the xyz file hardcoded into the subroutine
    !
    implicit none
    !
    !IO_variables
    !
    integer*2, intent(in):: natm
    real*8, intent(out) :: coords(3,natm)
    integer*1, intent(out) :: atmnum(natm)
    !
    ! Local variables
    !
    character(len=*), parameter :: file_path = "../../test_files/test_3b_1.xyz"
    integer :: ifi, ioflag, i
    integer*2 :: check_natm
    character :: symbol
    
    open(unit=ifi, file=file_path, err=9999, iostat=ioflag)
    if (ioflag .ne. 0) goto 9999

    read(ifi, *) check_natm
    if (check_natm .ne. natm) then
        print *, "Error, wrong number of files in xyz file. Was the right", &
                 "file used?"
        print*, "natoms should be ", natm
        stop 'read_xyzfile 2'
    endif

    read(ifi, *)
    do i=1, natm
        read(ifi, *, iostat=ioflag) symbol, coords(:,i)
        if (ioflag .ne. 0) then
            print *, "Error reading atom ", i
        endif
        if (symbol .eq. 'O') then
            atmnum(i) = 8
        elseif (symbol .eq. 'H') then
            atmnum(i) = 1
        else
            print *, "Could not identify atom"
            stop 'read_xyzfile 3'
        endif
    enddo

    return
9999 print *, "Error opening test file from path ", file_path
    stop 'read_xyzfile 1'
end subroutine read_xyzfile

