!> This a program that is used to call the symmetry function subroutines as
!! an external library. It uses the coordinates from one of the other test
!! runs
program caller
    use bp_symfuncs
    implicit none
    !
    ! Setup Variables
    !
    integer*1, parameter :: natm = 9
    real*8 :: coords(3,natm)
    integer*2 :: atmnum(natm)

    ! Read the .xyz coord file
    call read_xyzfile(natm, coords, atmnum)
    call print_coords(natm, coords, atmnum)

    ! This needs to always be called first to initialize the module
    call initialize_element_pars()
    ! Now we can define variables based on how it was initialized
    call stage_calculation(natm, coords, atmnum)


end program caller
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine stage_calculation(natm, coords, atmnum)
    ! After initialize_element_pars has been called, we can
    ! go and do the calculation
    use bp_symfuncs
    implicit none
    !
    ! IO Variables
    !
    integer*1 :: natm
    real*8 :: coords(3,natm)
    integer*2 :: atmnum(natm)
    !***********************************************
    ! BP Variables
    !***********************************************
    ! The basis functions to fill
    type(basis) rad_bas, ang_bas
    ! The number of geometries for the calculation
    integer, parameter :: num_geoms = 1
    ! The maximum number of atoms in the calculation
    integer :: max_atoms
    ! The number of each element. Think they are defined
    ! in order as H,O
    integer, parameter :: num_of_els(num_els) = [6,3]
    !
    type(mol_id), intent(inout) :: mol_ids(num_els)
    ! Below is index of molecules to basis functions. This is only
    ! Critical when multiple geometries are used
    integer*4, intent(inout) :: mol2bas(2, num_els, num_geoms)
    !***********************************************
    ! Local variables
    !

    max_atoms = int(natm)
    print *, 'hello'
    return

end subroutine stage_calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine print_coords(natm, coords, atmnum)
    ! Print the coordinates to stdout in xyz file format
    implicit none
    integer*1 :: natm
    real*8 :: coords(3,natm)
    integer*2 :: atmnum(natm)

    integer :: i

    print *, "Input coordinates:"
    print *, natm
    print *, 
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
    integer*1, intent(in):: natm
    real*8, intent(out) :: coords(3,natm)
    integer*2, intent(out) :: atmnum(natm)
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

