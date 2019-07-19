module libtfserver
    use iso_c_binding
    public :: tfserver

    ! Include the fortran-c interface function definitions
    !include "tfserving_utils_cdef.f90"
    ! C function definition for the TFServing object
    ! Functions are in tfserving_utils_capi.cpp
    interface
        ! Constructor
        function create_tfserver_c(port, model_name) bind(C, name="create_tfserver")
            use iso_c_binding
            type(c_ptr) :: create_tfserver_c
            integer(c_int), value :: port
            character(kind=c_char,len=1) :: model_name
        end function
        ! Destructor
        subroutine delete_tfserver_c(tfserver) bind (C, name="delete_tfserver")
            use iso_c_binding
            type(c_ptr), value :: tfserver
        end subroutine delete_tfserver_c
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

    subroutine delete_tfserver(this)
        implicit none
        type(tfserver) :: this
        call delete_tfserver_c(this%ptr)
    end subroutine

    ! Bounds procedure needs to take a polymorphic (class) argument
    subroutine delete_tfserver_polymorph(this)
        implicit none
        class(tfserver) :: this
        call delete_tfserver_c(this%ptr)
    end subroutine delete_tfserver_polymorph

end module libtfserver

