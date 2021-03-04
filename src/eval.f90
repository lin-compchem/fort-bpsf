! This program can be used to read a h5 coordinate file and talk to a
! tensorflow server server
!
#include "parameters.h"
program main                                                      
    use bp_symfuncs
    use libtfserver    
    use iso_c_binding
    use iso_fortran_env    
    use hdf5                                                                     
    use h5d                                                                      
    use h5_file_info
    !
    ! Setup Variables
    !
    implicit none                                                                
 
    character(len=:), allocatable :: if_path, of_path, h5_path, model_name                
    real*8, allocatable :: coords(:,:,:)                                         
    integer*ANUMKIND, allocatable :: atmnm(:,:)                                  
    integer*NATMKIND, allocatable :: natoms(:)
    type(tfserver), allocatable :: tfs
    !
    ! These are the variables we need for the model server
    !   
    integer(c_int32_t) :: port
    !
    ! Store the output of the model server
    !
    double precision, allocatable :: energy(:), gradient(:,:,:)
    ! Loop variables, etc.
    logical :: convert_to_au = .False.
    integer :: i, a
    integer*NATMKIND :: aa                                  
    double precision, parameter :: ang2au = 1.889725989
    double precision, parameter :: au2ang = 0.529177249
    double precision, parameter :: au2kcal = 627.5095d0
    double precision, parameter :: kcal2au = 1.00 / au2kcal
    double precision, parameter :: kcalpang2au = kcal2au/ang2au
    ! Get the paths from the command line                                        
    call get_args(if_path, of_path, port, model_name)
    call read_input_file(if_path, h5_path=h5_path)                               
    ! Get the data from the h5 file                                              
    call FI_init_from_ifi(h5_path)                                               
    call FI_get_data(coords, natoms, atmnm)
    call FI_close_input_h5()

    allocate(energy(num_geoms))
    allocate(gradient, mold=coords)

    allocate(tfs)
    tfs = tfserver(port, model_name)

    energy = 0d0
    gradient = 0d0
    
    do i=1, num_geoms
        a = natoms(i)
        aa = natoms(i)
        call tfs%wrapBPSF(coords(1:3,:a,i), aa, atmnm(:a,i), &
                          energy(i), gradient(1:3,:a,i))  
    enddo

    if (convert_to_au .eqv. .true.) then
        energy = energy * kcal2au
        gradient = gradient * kcalpang2au
        print *, "energy:  ", energy
        print *, "gradient"
        print *, gradient
    endif

    !
    ! Write the resulting file!
    !
    call write_engrad(energy, gradient, num_geoms, max_atoms, of_path)  
end program main 

subroutine write_engrad(energy, gradient, num_geoms, max_atoms, of_path)
    use h5d
    use hdf5
    use h5s
    double precision:: energy(num_geoms)
    double precision:: gradient(3,max_atoms,num_geoms)
    integer :: num_geoms, max_atoms
    integer(kind=4) :: error
    character(len=*), intent(in) :: of_path
    INTEGER, PARAMETER :: real_kind_15 = SELECTED_REAL_KIND(15,307)
    !
    ! HDF5 Variables and local
    !
    integer(HID_T) :: ofi ! output file handle
    integer(HID_T) :: dspace_id, dset_id   
    
    call h5fcreate_f(of_path, H5F_ACC_TRUNC_F, ofi, error)
    if (error .ne. 0)  goto 2000
    !
    ! write the energies
    !
    call h5screate_simple_f(1, [int8(num_geoms)], dspace_id, error)
    ! Create dataset
    CALL h5dcreate_f(ofi, 'energy', H5T_NATIVE_DOUBLE, &
                     dspace_id, dset_id, error)
    if (error .ne. 0) goto 1005
    ! Write data
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, energy(:), [int8(num_geoms)], error)
    if (error .ne. 0) goto 1010
    ! Close the dataset
    CALL h5dclose_f(dset_id , error)
    CALL h5sclose_f(dspace_id, error)
    if (error .ne. 0) goto 1015
    !
    ! write the gradients
    !
    call h5screate_simple_f(3, int8(shape(gradient)), dspace_id, error)
    ! Create dataset
    CALL h5dcreate_f(ofi, 'gradient', H5T_NATIVE_DOUBLE , dspace_id, dset_id, error)
    if (error .ne. 0) goto 1020
    ! Write data
    CALL h5dwrite_f(dset_id, H5T_NATIVE_double, gradient, int8(shape(gradient)), error)
    if (error .ne. 0) goto 1025
    ! Close the dataset
    CALL h5dclose_f(dset_id , error)
    CALL h5sclose_f(dspace_id, error)
    if (error .ne. 0) goto 1030
    !
    ! Close the file  
    !
    call h5fclose_f(ofi, error)
    return
    ! Error handling
1000 print *, "Error writing energy"
    stop 'save_basis 1'
1005 print *, "Error writing energy"
    stop 'save_basis 2'
1010 print *, "Error writing energy"
    stop 'save_basis 3'
1015 print *, "Error writing energy"
    stop 'save_basis 4'
1020 print *, "Error writing gradient"
    stop 'save_basis 5'
1025 print *, "Error writing gradient"
    stop 'save_basis 6'
1030 print *, "Error writing gradient"
    stop 'save_basis 7'
1035 print *, "Error writing energy"
    stop 'save_basis 8'
2000 print *, "error opening output file: ", of_path
    stop 'save_basis 17'
end subroutine write_engrad
