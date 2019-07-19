program fdemo1
    use libtfserver
    use iso_c_binding
    implicit none
    integer(c_int) :: port = 8504
    type(tfserver) :: tfs
    character(len=*), parameter :: model_name = "full_basis"
    ! create a new connection to the tensorflow server
    tfs = tfserver(port, model_name)
    ! GFortran doesn't always call destructor upon exit of scope
    !
    ! The below bug is fixed in gfortran >= 7.0 but I include it here because
    ! I don't know how to check gfortran version
    !
#ifdef __GNUC__
    call tfs%delete
#endif
end program fdemo1
