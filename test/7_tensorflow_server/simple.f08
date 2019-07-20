!> This a program that is used to call the symmetry function subroutines as
!! an external library. It uses the coordinates from one of the other test
!! runs
program caller
    use libtfserver
    use iso_c_binding
    implicit none
    integer(c_int) :: port = 8504 
    character(len=*), parameter :: model_name = "full_basis" 
    type(tfserver) :: tfs

    ! Initialize the TFServer Interface
    tfs = tfserver(port, model_name)
    
    ! Run the Test
    call tfs%modelTest1()

    call tfs%delete()
end program caller
