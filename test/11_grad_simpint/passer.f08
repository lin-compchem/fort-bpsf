!> This a program that is used to call the symmetry function subroutines as
!! an external library. It uses the coordinates from one of the other test
!! runs
#include "../../src/parameters.h"
program caller
    use bp_symfuncs
    use libtfserver
    use iso_c_binding
    implicit none
    !
    ! Setup Variables
    ! These in general require some user input but here they are
    ! hard coded
    !
    integer*NATMKIND, parameter :: natm = 9 ! for reading xyz file
    real*8 :: coords(3,natm)                ! coords from xyz file
    integer*ANUMKIND :: atmnum(natm)        ! atomic numbers
    integer(c_int) :: port = 8504           ! port from tfs script in test dir
    character(len=*), parameter :: model_name = "gradientmodel" ! name of model
                                            ! also from tfs script
    character(len=:), allocatable :: in_path! name of BPSF input file
    !
    ! These are variables for the interface that don't require input
    !
    type(tfserver) :: tfs ! object that contains interface
    real(c_double) :: energy ! This is the result of the calculation
    real(c_double) :: gradient(3,natm) ! Will be calculated in subroutine

    integer i
    in_path = 'bp.inp'
    energy = -1.0d0

    ! Check the TF server before things get going
    tfs = tfserver(port, model_name)

    ! Read the .xyz coord file
    call read_xyzfile(natm, coords, atmnum)
    call print_coords(natm, coords, atmnum)

    ! This needs to always be called first to initialize the BPSF module
    call read_input_file(in_path)

    ! Do the calculation
    call tfs%wrapBPSF(coords, natm, atmnum, energy, gradient)

    ! Print the result of the calculation
    print *, "The energy is: ", energy
    print *, "The energy should be: 11.927280400000001"
    
    if (energy -  11.927280400000001 .gt. 1e-6) then
        print *, "ERROR TEST FAILED"
    else
        print *, "TEST SUCCEEDED"
    endif
    
    ! print the gradient of the basis
    10 format( 3(F20.12,8x))
    do i=1, natm
        print 10, gradient(:,i)
    enddo 
    ! This step is optional for GCC >=7 which uses final procedure
    call tfs%delete()
end program caller
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine print_coords(natm, coords, atmnum)
    ! Print the coordinates to stdout in xyz file format
    implicit none
    integer*NATMKIND :: natm
    real*8 :: coords(3,natm)
    integer*ANUMKIND :: atmnum(natm)

    integer :: i, j

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
    integer*NATMKIND, intent(in):: natm
    real*8, intent(out) :: coords(3,natm)
    integer*ANUMKIND, intent(out) :: atmnum(natm)
    !
    ! Local variables
    !
    character(len=*), parameter :: file_path = "../../test_files/test_3b_1.xyz"
    integer :: ifi, ioflag, i
    integer*NATMKIND :: check_natm
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
