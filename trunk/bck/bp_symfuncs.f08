module bp_symfuncs
    implicit none
    integer, parameter :: dp = selected_real_kind(15,300)
    integer, parameter :: num_els = 2
    integer, parameter :: els(2) = [1, 8]

    ! General Variables
    real(dp) :: Rc = 6.5

    ! Radial variables
    integer, parameter :: max_etas = 50
    integer, parameter :: max_bonds = 2
    real(dp) :: eta(max_etas, max_bonds, num_els) = 0d0
    real(dp) :: rs(max_etas, max_bonds, num_els) = 0d0
    integer :: rad_size(max_bonds, num_els) = 0, radbas_length(num_els) = 0, &
     rad_types(max_bonds, num_els) = 0, num_bonds(num_els) = 0

    ! Angular Variables
    integer, parameter :: max_angles = 3
    integer, parameter :: max_lambdas = 2
    integer, parameter :: max_eprimes = 3
    integer, parameter :: max_zetas = 3
    real(dp) :: zetas(max_zetas, max_angles, num_els) = 0d0
    real(dp) :: eprimes(max_eprimes, max_angles, num_els) = 0d0
    real(dp) :: lambdas(max_lambdas, max_angles, num_els) = 0d0
    integer :: ang_size(max_angles, num_els) =0d0, angbas_length(num_els) = 0d0, &
               ang_types(2, max_angles, num_els),  num_angles(num_els) = 0

    ! integer, allocatable :: elements(:)
    ! real(dp), allocatable :: eprimes(:, :)
    ! integer, allocatable :: num_eta_rs(:), num_eprimes(:)
    ! integer, allocatable :: num_els
    ! integer, allocatable :: bonds(:)
    ! integer :: num_bonds
    ! integer :
    ! integer, allocatable :: angles(:, :)
    ! integer :: num_angles

contains
    subroutine initialize_rs_eta(rs_start, rs_end, num, rs, etas)
!    This subroutine (re)initializes the rs and eta arrays according to the
!    expression:
!        1 / sqrt(2*eta) = 0.2 * Rs
!    equal spacing is created between the rs_start and rs_end parameters
        implicit none
!        IO VARS
        real(dp), intent(in) :: rs_start, rs_end
        integer, intent(in) :: num
        real(dp), intent(inout):: rs(:), etas(:)
!        Local Vars
        real(dp) rs_spacer
        integer i
!        Begin
        rs_spacer = (rs_end - rs_start) / (num - 1)
        do i=0, num - 1
            rs(i + 1) = rs_start + i * rs_spacer
        enddo
        etas(1:num) = 0.5 * ( 5. / Rs(1:num)) ** 2
    end subroutine initialize_rs_eta
    
    subroutine initialize_element_pars
        ! Initialize the element parameters. This can be modified to add different ones in
        implicit none
        integer i, j, k, l, m
        real(dp) :: rs_start = 0.8d0, rs_end = 8.0d0
        integer :: num_rs = 24
        ! Initialize all vars with 0s

        
        ! Initialize the eta RS values
        rad_types(1, 1) = 1
        rad_size(1, 1) = num_rs
        call initialize_rs_eta(rs_start, rs_end, num_rs, rs(:, 1, 1), eta(:, 1, 1))
        rad_types(2, 1) = 8
        rad_size(2, 1) = num_rs
        call initialize_rs_eta(rs_start, rs_end, num_rs, rs(:, 2, 1), eta(:, 2, 1))
        num_bonds(1) = 2
        
        rad_types(1, 2) = 1
        rad_size(1, 2) = num_rs
        call initialize_rs_eta(rs_start, rs_end, num_rs, rs(:, 1, 2), eta(:, 1, 2))
        rad_types(2, 2) = 8
        rad_size(2, 2) = num_rs - 2
        call initialize_rs_eta(1.4260869565217393d0, rs_end, 22, rs(:, 2, 2), eta(:, 2, 2))
        num_bonds(2) = 2

        do i=1, num_els
            radbas_length(i) = 0d0
            do j=1, max_bonds
                radbas_length(i) = radbas_length(i) +rad_size(j, i)
            enddo
        enddo

        !Initialize the angle variables
        ang_types(:, 1, 1) = [1, 1]
        ang_types(:, 2, 1) = [1, 8]
        num_angles(1) = 2
        ang_types(:, 1, 2) = [8, 8]
        ang_types(:, 2, 2) = [1, 8]
        ang_types(:, 3, 2) = [1, 1]
        num_angles(j) = 3

        do i=1, num_els
            angbas_length(i) = 0
            do j=1, num_angles(i)
                lambdas(1:2, j, i) = [-1.d0, 1.d0]
                eprimes(1:3, j, i) = [0.001d0, 0.01d0, 0.05d0]
                zetas(1:3, j, i) = [1.d0, 4.d0, 16.d0]
                ang_size(j, i) = 2 * 3 * 3
                angbas_length(i) = angbas_length(i) + ang_size(j, i)
            enddo
        enddo
        
        
    end subroutine initialize_element_pars
end module bp_symfuncs