#include "../parameters.h"
module libtfserver
    use iso_c_binding
    implicit none
    public :: tfserver

    ! Include the fortran-c interface function definitions
    !include "tfserving_utils_cdef.f90"
    ! C function definition for the TFServing object
    ! Functions are in tfserving_utils_capi.cpp
    interface
        !
        ! Constructor
        function create_tfserver_c(port, model_name) bind(C, name="create_tfserver")
            use iso_c_binding
            type(c_ptr) :: create_tfserver_c
            integer(c_int), value :: port
            character(kind=c_char,len=1) :: model_name
        end function
        !
        ! Destructor
        subroutine delete_tfserver_c(tfserver) bind (C, name="delete_tfserver")
            use iso_c_binding
            type(c_ptr), value :: tfserver
        end subroutine delete_tfserver_c
        !
        ! Run a test with a simplified TF model
        subroutine model_test1_c(tfserver) bind(C, name="tfs_model_test1")
            use iso_c_binding
            type(c_ptr), value :: tfserver
        end subroutine model_test1_c
        ! Print the basis
        subroutine bpsf_energy_c(tfserver, bas, max_bas, max_atoms, &
            num_bas, num_of_els, num_els, energy) bind(C, name="tfs_bpsf_energy")
            use iso_c_binding
            type(c_ptr), value :: tfserver
            real(c_double), intent(in) :: bas(max_bas, max_atoms, num_els)
            integer(c_int), intent(in) :: max_bas         ! dimension of basis array
            integer(c_int), intent(in) :: max_atoms       ! dimension of basis array
            integer(c_int), intent(in) :: num_bas(num_els)! number of basis functions for each el
            integer(c_int), intent(in) :: num_els         ! total number elements
            integer(c_int), intent(in) :: num_of_els(num_els)! number of each element
            real(c_double), intent(out) :: energy ! energy of calculation
        end subroutine bpsf_energy_c
    end interface

    ! Use a fortran type to represent the c++ class in an opaque manner
    type tfserver
        type(c_ptr) :: ptr ! pointer to the TFServer class
    contains
    ! "We can bind some functions to this type allowing for a cleaner
    !  syntax"
    ! Required because
    ! GFortran doesn't always call destructor upon exit of scope
    !
    ! The below bug is fixed in gfortran >= 7.0 but I include it here because
    ! I don't know how to check gfortran version
    !
#ifdef __GNUC__
        procedure :: delete => delete_tfserver_polymorph
#else
        final :: delete_tfserver ! Destructor
#endif
        procedure :: sendBPSF => bpsf_energy
        procedure :: modelTest1 => model_test1
    end type tfserver

    ! This function acts as a constructor for the tfserver type
    interface tfserver
        procedure create_tfserver
    end interface

contains ! Implementation of the functions, we just wrap the C function here
    function create_tfserver(port, model_name) ! Constructor
        implicit none
        type(tfserver) :: create_tfserver
        integer, intent(in) :: port
        character(len=*), intent(in) :: model_name
        create_tfserver%ptr = create_tfserver_c(port, &
                                    trim(adjustl(model_name))//char(0))
    end function create_tfserver

    subroutine delete_tfserver(tfs)
        implicit none
        type(tfserver) :: tfs
        call delete_tfserver_c(tfs%ptr)
    end subroutine

    ! Bounds procedure needs to take a polymorphic (class) argument
    subroutine delete_tfserver_polymorph(tfs)
        implicit none
        class(tfserver) :: tfs
        call delete_tfserver_c(tfs%ptr)
    end subroutine delete_tfserver_polymorph

    ! subroutine bpsf_energy(this)
    subroutine bpsf_energy(tfs, basis, max_bas, max_atoms, &
            num_bas, num_of_els, num_els, energy)
        class(tfserver), intent(in) :: tfs
        real(kind=8), intent(in) :: basis(max_bas, max_atoms, num_els)
        integer(kind=4), intent(in) :: max_bas         ! dimension of basis array
        integer(kind=4), intent(in) :: max_atoms       ! dimension of basis array
        integer(kind=4), intent(in) :: num_bas(num_els)! number of basis functions for each el
        integer(kind=4), intent(in) :: num_of_els(num_els)! number of each element
        integer(kind=4), intent(in) :: num_els         ! total number elements
        real(kind=8), intent(out) :: energy ! energy of calculation

        call bpsf_energy_c(tfs%ptr, basis, max_bas, max_atoms, &
           num_bas, num_of_els, num_els, energy)
    end subroutine bpsf_energy

    subroutine model_test1(tfs)
        implicit none
        class(tfserver) :: tfs
        call model_test1_c(tfs%ptr)
    end subroutine model_test1

    ! This subroutine is a wrapper subroutine that accepts coordinates and
    !   returns whatever
    subroutine wrap_bpsf(tfs, coordinates, num_atoms, atomic_numbers, energy, &
            gradient)
        use bp_symfuncs, only: calc_bp, num_els, max_bas, radbas_length, angbas_length
        !
        ! I/O Variables
        !
        class(tfserver), intent(in) :: tfs
        integer*NATMKIND, intent(in) :: num_atoms
        integer*ANUMKIND, intent(in) :: atomic_numbers(num_atoms)
        real(kind=8), intent(in) :: coordinates(3, num_atoms)
        real(kind=8), intent(out) :: energy
        real(kind=8), intent(out), optional :: gradient
        !
        ! Local Variables
        !
        integer :: g_num_of_els(2) = [0,0]
        real*8, dimension(max_bas, num_atoms, num_els) :: rad_bas, ang_bas
        real*8, dimension(3, num_atoms, max_bas, num_atoms, num_els) :: rad_grad, ang_grad
        integer el
        integer(kind=4) max_atoms
        max_atoms = num_atoms
        !integer el
        call calc_bp(num_atoms, coordinates, atomic_numbers, rad_bas, ang_bas, &
                     rad_grad, ang_grad, max_atoms, g_num_of_els)        

        ! Concatenate the basis functions into the radbas array
        do  el=1, num_els
           rad_bas(radbas_length(el)+1:radbas_length(el)+angbas_length(el),:g_num_of_els(el),el) = &
               ang_bas(:angbas_length(el),:g_num_of_els(el),el)
        enddo
    
        ! Send the concated basis to the server
        call tfs%sendBPSF(rad_bas, max_bas, max_atoms, radbas_length(:) + &
            angbas_length(:), g_num_of_els, num_els, energy)
  
        if (present(gradient)) then
            stop "I shouldn't be here"
        endif
    end subroutine wrap_bpsf
end module libtfserver

