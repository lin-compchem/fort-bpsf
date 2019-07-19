! C function definition for the TFServing object
! Functions are in tfserving_utils_capi.cpp
use iso_c_binding
interface
    ! Constructor
    function create_tfserver_c(port) bind(C, name="create_tfserver")
        type(c_ptr) :: create_tfserver_c
        integer(c_int), value :: port
    end function
    ! Destructor
    subroutine delete_tfserver_c(tfserver) bind (C, name=delete_tfserver)
        type(c_ptr), value :: tfserver
end interface
