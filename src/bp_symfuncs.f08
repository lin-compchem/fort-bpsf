module bp_symfuncs
    implicit none
    integer, parameter :: dp = selected_real_kind(15,300)

    real(8), parameter :: pi = 4 * atan(1.0d0)
    integer, parameter :: max_bas = 150
    logical, parameter :: calc_grad = .true.
    integer, parameter :: debug = 1
    !TODO: Check for good value for max_bas

    ! General Variables
    real(dp) :: Rc = 0.0                  ! Rc cutoff in bp functions
    integer :: num_els = 0                ! Number of elements to calc bf for
    integer, parameter :: max_bonds = 0   ! Maximum number of BOND TYPES (categories of rad functions)
    integer, parameter :: max_etas = 0    ! Maximum pairs of rs_eta
    integer, parameter :: max_angles = 0  ! Maximum number of ANGLE TYPES (categories of ang functions)
    integer :: max_ezl = 0              ! Maximum triplets of eprime/zeta/lambda


    real(dp):: pi_div_Rc = 0.0d0
    integer, allocatable :: els

    ! Radial variables
    real(dp), allocatable :: eta(:,:,:)
    real(dp), allocatable :: rs(:,:,:)

    integer, allocatable :: &
       rad_size(:,:), & ! The size of each radial basis set (max_eta for each atom type)
       radbas_length(:), & ! The total lenght of the radial basis for each element
       rad_types(:,:), & ! The radial basis bond types (atomic numbers)
       num_bonds(:), &  ! The number of bond types for each element
       rad_b_ind(:,:) ! ? Think this is where in the radbas each bond starts

    ! Angular Variables
    real(dp), allocatable :: &
        etzetlam(:,:,:,:)  ! The eta/zeta/lambda variables. Array goes ([e,z,l],max_el,ang_type,el)
    integer, allocatable :: &
        ang_size(:,:), & ! The number of angular basis functions to calculate for each angle type
        angbas_length(:), & ! The total angular basis length for each element
        ang_types(:,:,:), & ! The ang types for each element
        num_angles(:),  & ! The number of ang types for each element
        ang_b_ind(:,:) ! ? Think this is where each angbas starts in the main basis


    ! Timing variables
    logical :: dotime = .false.
    integer :: timing_interval = 100
    
    integer(kind=8), parameter :: inf = 789555466321 !input file handle

    type basis
        ! The basis functions. First dimension corresponds to atom, second dim
        ! is the basis value
        real*8, dimension(:,:), allocatable :: b
        ! The gradient of the basis in BxExNx3 dimensions. The first is the
        ! number of basis funcitons for a given atom, the second has N dimensions
        ! where N is the maxbas parameter and holds the atoms for a. The thirs
        real*8, dimension(:,:,:,:), allocatable :: g
    end type basis

    type mol_id
        ! This links the atom basis to the molecule. There is one entry for
        ! each atom of a given element. 0 would correspond to belonging to the
        ! first molecule, 1 the second, etc.
        integer*4, dimension(:), allocatable :: bas2mol
    end type mol_id


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
        ! Initialize the element parameters (num bonds, etas, zetas, etc.)
        implicit none
        integer i, j
        real(dp) :: rs_start = 0.8d0, rs_end = 8.0d0
        integer :: num_rs = 24
        ! Initialize all vars with 0s

        
        ! Initialize the eta RS values
        rad_types(1, 1) = 1
        rad_size(1, 1) = num_rs
        rad_b_ind(1, 1) = 1
        call initialize_rs_eta(rs_start, rs_end, num_rs, rs(:, 1, 1), eta(:, 1, 1))
        rad_types(2, 1) = 8
        rad_size(2, 1) = num_rs
        rad_b_ind(2, 1) = rad_b_ind(1, 1) + rad_size(1, 1)
        call initialize_rs_eta(rs_start, rs_end, num_rs, rs(:, 2, 1), eta(:, 2, 1))
        num_bonds(1) = 2
        
        rad_types(1, 2) = 1
        rad_size(1, 2) = num_rs
        rad_b_ind(1, 2) = 1
        call initialize_rs_eta(rs_start, rs_end, num_rs, rs(:, 1, 2), eta(:, 1, 2))
        rad_types(2, 2) = 8
        rad_size(2, 2) = num_rs
        rad_b_ind(2, 2) = rad_b_ind(1, 2) + rad_size(1, 2)
        call initialize_rs_eta(rs_start, rs_end, num_rs, rs(:, 2, 2), eta(:, 2, 2))
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
        num_angles(2) = 3

        ang_b_ind(1, :) = 0
        do i=1, num_els
            angbas_length(i) = 0
            do j=1, num_angles(i)
                ang_b_ind(j, i) = ang_b_ind(j, i) + angbas_length(i) + 1
                lambdas(1:2, j, i) = [-1.d0, 1.d0]
                eprimes(1:3, j, i) = [0.001d0, 0.01d0, 0.05d0]
                zetas(1:3, j, i) = [1.d0, 4.d0, 16.d0]
                ang_size(j, i) = 2 * 3 * 3
                angbas_length(i) = angbas_length(i) + ang_size(j, i)
            enddo
        enddo
        !TODO: THIS MUST BE CHANGED IF THE ETAS ZETAS LAMBDAS ARE NOT ALL
        !      THE SAME
        call init_eta_zeta_lambda(eprimes(:,1,1), zetas(:,1,1), lambdas(:,1,1),&
             num_etzetlam, etzetlam)
    end subroutine initialize_element_pars

    subroutine init_eta_zeta_lambda(in_eta, in_zeta, in_lambda, &
                                    o_num_ezl, o_ezl)
    ! Initialize the ezl array which has all possible etas, zetas, and
    ! lambdas stored in an array for vectorized computations
    implicit none
    ! I/O Variables
    real(kind=8), intent(in) :: in_eta(:), in_zeta(:), in_lambda(:)
    integer, intent(out) :: o_num_ezl
    real(kind=8), allocatable, intent(inout) :: o_ezl(:,:)
    ! Local variables
    integer :: my_neta, my_nzeta, my_nlambda
    integer :: e, z, l, i

    if (rank(in_eta) .gt. 1) then
        stop 'init_eta_zeta_lambda invalid eta rank'
    elseif (rank(in_zeta) .gt. 1) then
        stop 'init_eta_zeta_lambda invalid zeta rank'
    elseif (rank(in_lambda) .gt. 1) then
        stop 'init_eta_zeta_lambda invalid lambda rank'
    endif

    my_neta = size(in_eta)
    my_nzeta = size(in_zeta)
    my_nlambda = size(in_lambda)
    o_num_ezl = my_neta * my_nzeta * my_nlambda
    allocate(o_ezl(o_num_ezl, 3))
    i = 1
    do e=1, num_eprimes
        do z=1, num_zetas
            do l=1, num_lambdas
                o_ezl(i, 1:3) = [in_eta(e), in_zeta(z), in_lambda(l)]
                i = i + 1
            ENDDO
        ENDDO
    ENDDO

    end subroutine init_eta_zeta_lambda

    subroutine calculate_basis(rad_bas, ang_bas, coords, atmnms, natoms, &
                               max_atoms, num_geoms, num_of_els, mol_ids, &
                               mol2bas)
        ! The meat of the program, lets calculate the basis functions!
        implicit none
        ! I/O vars
        integer, intent(in) :: num_geoms, max_atoms
        integer, intent(in) :: num_of_els(num_els)
        real*8, intent(in) :: coords(3, max_atoms, num_geoms)
        integer*1, intent(in) :: atmnms(max_atoms, num_geoms)
        integer*2, intent(in) :: natoms(num_geoms)
!TODO:  MAKE SURE THIS CORRESPONDS WITH THE ACTUAL DIMS IF TOO MUCH MONKEYING IS DONE
        type(basis), intent(inout) :: rad_bas(num_els), ang_bas(num_els)
        type(mol_id), intent(inout) :: mol_ids(num_els)
        integer*4, intent(inout) :: mol2bas(2, num_els, num_geoms)
        !
        ! Local vars
        !
        ! assorted counters
        integer :: i, g
        ! g : the current geometry


        !Basis sets to be reduced in openmp loops
        real*8, dimension(max_bas, max_atoms, num_els) :: tmp_rad_bas, tmp_ang_bas
        real*8, dimension(3, max_atoms, max_bas, max_atoms, num_els) :: tmp_rad_grad, tmp_ang_grad
        ! Number of atoms for a given element for a geometry
        integer :: g_num_of_els(num_els)
        integer :: natm
        ! Counters to keep track of elements location in basis set throughout
        ! multiple geometries.
        integer :: bas_s(num_els), i_end
        ! Variables for timing
        real(kind=8) :: t
        integer*8 :: begin_time, begin_loop

        bas_s(:) = 1
        !!!! Con
        rad_bas(1)%g(:,:,:,:) = 0.d0
        rad_bas(2)%g(:,:,:,:) = 0.d0
        ang_bas(1)%b(:,:) = 0.d0
        ang_bas(2)%b(:,:) = 0.d0
        ang_bas(1)%g(:,:,:,:) = 0.d0
        ang_bas(2)%g(:,:,:,:) = 0.d0

        ! Calculate the vectors and rijs
        call tick(begin_time)
        call tick(begin_loop)
    over_geoms: do g=1, num_geoms
        natm = natoms(g)
        call calc_bp(natoms(g), coords(:,:natm,g), atmnms(:natm,g), &
                     tmp_rad_bas, tmp_ang_bas, &
                     tmp_rad_grad,tmp_ang_grad, &
                     max_atoms, g_num_of_els)

! Copy the temporary basis sets to our main basis
save_basis: do i=1, num_els
            i_end = bas_s(i) + g_num_of_els(i) - 1
            rad_bas(i)%b(:,bas_s(i):i_end) = tmp_rad_bas(:radbas_length(i), &
                :g_num_of_els(i), i)
            rad_bas(i)%g(:,:,:,bas_s(i):i_end) = &
                tmp_rad_grad(:, :, :radbas_length(i), :g_num_of_els(i), i)
            ang_bas(i)%b(:,bas_s(i):i_end) = tmp_ang_bas(:angbas_length(i), &
                :g_num_of_els(i), i)
            ang_bas(i)%g(:,:,:,bas_s(i):i_end) = &
                tmp_ang_grad(:, :, :angbas_length(i), :g_num_of_els(i), i)
            mol_ids(i)%bas2mol(bas_s(i):i_end) = g - 1
            mol2bas(1:2,i,g) = [bas_s(i)-1, i_end]
            bas_s(i) = i_end + 1
        enddo save_basis


     ! Its time for time
        if (mod(g, timing_interval) .eq. 0) then
            t = tock(begin_loop)
            print *, 'Loop iterations', timing_interval, 'Time: ', t
            call tick(begin_loop)
        endif
    enddo over_geoms
    t = tock(begin_time)
    print *, "Total time: ", t

    return
    end subroutine calculate_basis

pure subroutine get_triplets(max_ind, max_trip, triplets, num_triplets)
!
! NOT WORKING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DONT USE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Get the number of unique triplets from 1 to max_ind. Return the triplets in
! the triplets variable and return the number of triplets in the num_triplets
! variable
!
    implicit none
    ! I/O variables
    integer, intent(in) :: max_ind
    integer, intent(in) :: max_trip
    integer, intent(out) :: triplets(3, max_trip)
    integer, intent(out) :: num_triplets
    ! local variables
    integer :: i, j, k

    num_triplets = 0
    do i=1, max_ind
        do j = 1, max_ind - 1
            if ( i .ne. j) then
                do k=j+1, max_ind
                    if ( i .eq. k) then
                        continue
                    elseif (j .eq. k) then
                        continue
                    else
                        num_triplets = num_triplets + 1
                        triplets(:,num_triplets) = [i, j, k]
                    endif
                enddo
            else
                continue
            endif
        enddo
    enddo
end subroutine get_triplets

pure subroutine calc_ij_for_ut(i, j, x, natm)
! Calculate the i, j coordiate for a square matrix based upon the
! 0-baSED index of the upper triangle. For example, for a three by three matrix,
! the indicies would range from 0-2 corresponding to the ij pairs (0,1);(0,2);(1,2)
! This will return the ij pairs (1,2);(1,3);(2,3)
    integer, intent(out) :: i, j
    integer, intent(in) :: x, natm
    i = x/natm
    j = mod(x, natm)
    if (j .le. i) then
        i = natm - i - 2
        j = natm - j - 1
    endif
    i = i + 1
    j = j + 1
end subroutine

subroutine tick(t)
    integer*8, intent(OUT) :: t

    call system_clock(t)
end subroutine tick

! returns time in seconds from now to time described by t
real(kind=dp) function tock(t)
    integer*8, intent(in) :: t
    integer*8 :: now, clock_rate

    call system_clock(now,clock_rate)

    tock = real(now - t)/real(clock_rate)
end function tock

subroutine find_num_els(natm, atmnms, el_key, g_num_of_els, &
    atom_basis_index)
! Go throgh the atoms in the geometry and find out what element each
! atom is. Then tally the number of each element.
! We output two things in this subroutine.
!    1. el_key:
!            The total number of each element.
!    2. g_num_of_els:
!           This subroutine only adds to this so it must be initialized before
!           this routine
!    3. atom_basis_index
!           An array with len natm that is not the cumulative number of elements
!           I.e. the vector of elements [O H H H O H H] would be:
!                                        1 1 2 3 2 4 5
    implicit none
    !
    ! IO Vars
    !
    integer, intent(in) :: natm
    integer*1, intent(in) :: atmnms(natm)
    integer, intent(inout) :: el_key(natm), & ! The element type that each atom belongs to
                   g_num_of_els(num_els), & ! The number of atoms for each element
                   atom_basis_index(natm)
    !
    ! Local Vars
    !
    integer i, j
    logical el_found
    !
    ! Reinitialize variables
    !
    g_num_of_els = 0
 do_natoms: do i =1, natm
                el_found = .false.
                do j=1, num_els
                    if (atmnms(i) .eq. els(j)) then
                        g_num_of_els(j) = g_num_of_els(j) + 1
                        el_found = .true.
                        el_key(i) = j
                        atom_basis_index(i) = g_num_of_els(j)
                    endif
                enddo
                ! Error handling
                if (el_found .eqv. .false.) then
                    print *, 'Error, could not find element in geom for atom',i
                    stop 'calc_bas 1'
                endif
           enddo do_natoms
end subroutine find_num_els

! 'standalone' basis function calculation. Here it is abstracted
! from the h5ls stuff and contains the raw data needed to calculate
! the basis functions for 1 geometry.
subroutine calc_bp (natm, coords, atmnms, rad_bas, ang_bas, &
                    rad_grad, ang_grad, max_atom, g_num_of_els)
    implicit none
    !
    ! Input Variables
    !
    integer*2, intent(in) :: natm           ! number of atoms
    integer*1, intent(in) :: atmnms(natm) ! atomic numbers
    real*8, intent(in) :: coords(3, natm) ! coordinates
    integer, intent(in) :: max_atom ! size of basis
    !
    !Basis sets to be reduced in openmp loops
    !
    real*8, dimension(max_bas, max_atom, num_els), intent(inout) :: &
                                            rad_bas, ang_bas
    real*8, dimension(3, max_atom, max_bas, max_atom, num_els), intent(inout) :: &
                                            rad_grad, ang_grad
    !
    ! Local Variables:
    !
    real*8, dimension(natm, natm) :: rij
    real*8, dimension(3, natm, natm) :: drijdx
    integer ut      ! upper triangle counter
    integer max_cos ! counter for unrolled cos matrix
    integer num_cos ! vestigial integer for if max_cos >= natm
    integer :: cos_inds(3, natm**2 * (natm-1)/2 )
    ! These also correspond to number of triplets, maximum number of triplets
    !
    ! Various counters
    integer a, b, c, e
    integer i, j, k, l, m
    integer x
    ! atom types for the ith and jth atom
    integer i_btype, j_btype, my_type

    ! Smoothing function array. fc(i,j) = fc(j,i) uses cosine func.
    !  dfcdr == dfc(rij)/dRij
    real*8, dimension(natm, natm) :: fc
    real*8, dimension(3, natm, natm) :: dfcdx

    ! You can look up an atom index and you will get 1 for H, 2 for O
    integer :: el_key(natm)
    ! Number of atoms for a given element for a geometry
    integer :: g_num_of_els(num_els)

    ! The index in the final basis set for a given atom. I.e. the first oxygen
    ! analyzed will have index of 1, second index of 2, ... then the first hydrogen
    ! will have index of one, etc...
    integer :: atom_basis_index(natm)
    !Int*4 version of natom
    integer :: natom

    ! Holds the radial basis during vectorized computation.
    real*8 :: tmp_rad(max_etas), tmp_drad(3, max_etas, 2), term,&
            tmp_dang(3, num_etzetlam, 3)
    ! Intermediary variables for calculating angular basis functions
    real*8 :: mycos, myexp, myfc, myang(num_etzetlam), mydcos(3,3)
    real*8 :: mya(num_etzetlam), myb(num_etzetlam), &
              ang_coeff(num_etzetlam)
    real*8 :: mydfc(3,3), mydexp(3,3)
    ! Get the number of upper triangular indicies
    ut = natm*(natm-1)/2
    ! The maximum number of cosines to calculate is n upper triangles
    max_cos = natm * ut
    g_num_of_els(:) = 0
    atom_basis_index(:) = 0
    ! This should be the same between calculations
    ang_coeff(:) = 2 ** (1 - etzetlam(:,2))
    ! set the diagonals of some arrays to 0
    drijdx(:,:,:) = 0d0
    natom = int(natm,kind=4)
    rad_bas = 0.0d0
    rad_grad = 0.0d0
    ang_bas = 0.0d0
    ang_grad = 0.0d0

    call find_num_els(natom, atmnms, el_key, g_num_of_els, atom_basis_index)
         ! generate the cosine inds. This is an array of all of the unique
         ! angles in the system. This is the key for the threads to use
         ! when the parallel parts start.
    call get_triplets(natom, max_cos, cos_inds, num_cos)

   !$OMP PARALLEL &
   !$OMP& private(i_btype, j_btype, tmp_rad) &
   !$OMP& private(x, i, j, k, a, b, c, e, m) &
   !$OMP& private(my_type, mycos, myexp, myfc, myang) &
   !$OMP& private(term, tmp_drad, mydcos, mydfc, mydexp) &
   !$OMP& private(mya, myb, tmp_dang)
   ! Remove from reduction and just add to critical cause we're allready doing
   ! it...
   !!$OMP& reduction(+:rad_bas) reduction(+:ang_bas)
   ! Cannot reduce the below arrays because they are too big?
   !!$OMP& reduction(+:rad_grad) reduction(+:ang_grad)
   !$OMP    DO
           ! Calculate the r-ij vectors and distances and cutoff

calc_bonds: do x=0, ut - 1
            call calc_ij_for_ut(i, j, x, natom)

            ! Distances
            drijdx(:, i, j) =  coords(:, i) - coords(:, j)
            rij(i, j) = norm2(drijdx(:, i, j))
            drijdx(:, i, j) = drijdx(:, i, j) / rij(i, j)

            ! Calculate the other diagonal
            rij(j, i) = rij(i, j)
            drijdx(:, j, i) = -drijdx(:, i, j)

            ! Cutoff
            if (rij(i, j) .lt. Rc) then
                fc(i, j) = 0.5 * ( cos(pi_div_Rc * rij(i, j)) + 1)
                term = -0.5 * pi_div_Rc * sin(pi_div_RC * rij(i,j))
                dfcdx(:, j, i) = term * drijdx(:, j, i)
                dfcdx(:, i, j) = term * drijdx(:, i, j)
            else
                fc(i, j) = 0d0
                dfcdx(:, i, j) = 0d0
                dfcdx(:, j, i) = 0d0
            endif
            fc(j, i) = fc(i, j)
            ! Find the bond type for the pair
            i_btype = 0
            j_btype = 0
            do k=1, num_bonds(el_key(i))
                if (rad_types(k, el_key(i)) .eq. atmnms(j)) then
                    i_btype = k
                    exit
                endif
            enddo
            do k=1, num_bonds(el_key(j))
                if (rad_types(k, el_key(j)) .eq. atmnms(i)) then
                    j_btype = k
                    exit
                endif
            enddo

            !Below should be changed if the etas change significantly as it
            !Assumes the etas are the same for all elements to save time
            ! This cannot be reordered because we use intermediary calculation
            ! Results for the variables
            tmp_rad(:) = exp(-1 * eta(:,1,1) * (rij(i,j) - rs(:,1,1))**2)
            ! Remove this later
            myexp = tmp_rad(1)
            ! The derivative calculation
            do k=1, 2
             do l=1, 3
                tmp_drad(l,:,k) = -2 * eta(:,1,1) * (rij(i,j) - rs(:,1,1)) * fc(i, j)
             enddo
            enddo
            do k=1, max_etas
                tmp_drad(:,k,1) = tmp_drad(:,k,1) * drijdx(:, i, j) + dfcdx(:, i, j)
                tmp_drad(:,k,2) = tmp_drad(:,k,2) * drijdx(:, j, i) + dfcdx(:, j, i)
            enddo
            do k=1, 2
             do l=1, 3
                tmp_drad(l,:,k) = tmp_drad(l,:,k) * tmp_rad(:)
             enddo
            enddo

            tmp_rad(:) = tmp_rad(:) * fc(i,j)

            if (i_btype > 0) then
                a = rad_b_ind(i_btype, el_key(i))
                b = a + rad_size(i_btype, el_key(i)) - 1
                c = atom_basis_index(i)
                e = el_key(i)

                !$OMP CRITICAL
                rad_bas(a:b, c, e) =  rad_bas(a:b, c, e) + tmp_rad(:)
                rad_grad(:, i, a:b, c, e) = rad_grad(:, i, a:b, c, e) + tmp_drad(:, :, 1)
                rad_grad(:, j, a:b, c, e) = rad_grad(:, j, a:b, c, e) + tmp_drad(:, :, 2)
                !$OMP END CRITICAL
            endif
            !
            ! Here i am adding the mirror RBF and corresponding gradients. Note that
            ! The assumption is that there is a corresponding basis function between H and O
            ! and that the gradient of that basis function is antisymmetric with it's
            ! corresponding gradient.
            !
            if (j_btype > 0) then
                a = rad_b_ind(j_btype, el_key(j))
                b = a + rad_size(j_btype, el_key(j)) - 1
                c = atom_basis_index(j)
                e = el_key(j)
                !$OMP CRITICAL
                rad_bas(a:b, c, e) =  rad_bas(a:b, c, e) + tmp_rad(:)
                rad_grad(:, i, a:b, c, e) = rad_grad(:, i, a:b, c, e) - tmp_drad(:, :, 1)
                rad_grad(:, j, a:b, c, e) = rad_grad(:, j, a:b, c, e) - tmp_drad(:, :, 2)
                !$OMP END CRITICAL
            endif
            enddo calc_bonds
    !$OMP ENDDO
!
! Main loop over angles
!
    !$OMP DO
calc_ang: do x=1, num_cos
            i = cos_inds(1, x)
            j = cos_inds(2, x)
            k = cos_inds(3, x)

            ! Find which angle type this triplet corresponds to
            my_type = 0
            do m=1, num_angles(el_key(i))
                if (ang_types(1, m, el_key(i)) .eq. atmnms(j)) then
                    if (ang_types(2, m, el_key(i)) .eq. atmnms(k)) then
                        my_type = m
                    endif
                elseif (ang_types(2, m, el_key(i)) .eq. atmnms(j)) then
                    if (ang_types(1, m, el_key(i)) .eq. atmnms(k)) then
                        my_type = m
                    endif
                endif
            enddo
            if (my_type .eq. 0) cycle
            !
            ! Calculate the cosine term
            !
            mycos = (rij(i,j)**2 + rij(i,k)**2 - rij(j,k)**2) / &
                    (2.0 * rij(i,j) * rij(i,k))
            ! Derivative of cos with respect to j, k. 2 is j, 3 is k
            term = 1/(rij(i,j) * rij(i,k))
            ! Derivative of cos in the j direction
            mydcos(:,2) = -mycos / rij(i,j) * drijdx(:,j,i)
            mydcos(:,2) = mydcos(:,2) + term * (coords(:,k) - coords(:,i))
            ! Derivative of cos in the k direction
            mydcos(:,3) = -mycos / rij(i,k) * drijdx(:,k,i)
            mydcos(:,3) = mydcos(:,3) + term * (coords(:,j) - coords(:,i))
            ! Derivative of cos with respect to i
!            mydcos(:,1) = -mycos * (rij(i,k) * drijdx(:,i,j) + rij(i,j)*drijdx(:,i,k))
!            mydcos(:,1) = term * (mydcos(:,1) + 2*coords(:,i,g) - coords(:,j,g) - coords(:,k,g))
            ! Actually, the derivative is equal and opposite to the sum of dj and dk
            mydcos(:,1) = -mydcos(:,2)-mydcos(:,3)

            ! Derivatives of fc(ij)*fc(ik)*fc(jk)
            mydfc(:,1) = fc(j,k)*( fc(i,j)*dfcdx(:,i,k) + fc(i,k)*dfcdx(:,i,j) )
            mydfc(:,2) = fc(i,k)*( fc(i,j)*dfcdx(:,j,k) + fc(j,k)*dfcdx(:,j,i) )
            mydfc(:,3) = fc(i,j)*( fc(i,k)*dfcdx(:,k,j) + fc(j,k)*dfcdx(:,k,i) )

            ! the myexp term is written with 2 different expressions. this is
            ! from the 2018 paper
            myexp = (rij(i, j) + rij(i, k) + rij(j, k))**2
            myfc = fc(i, j) * fc(i, k) * fc(j, k)

            ! Begin
            ! The gaussian parts that don't depend on eta
            mydexp(:,1) = drijdx(:,i,j) + drijdx(:,i,k)
            mydexp(:,2) = drijdx(:,j,i) + drijdx(:,j,k)
            mydexp(:,3) = drijdx(:,k,i) + drijdx(:,k,j)

            !
            ! Myang is the calculation for the angular basis functions
            !
            mya(:) = ( 1 + etzetlam(:,3) * mycos) ** etzetlam(:,2)
            myb(:) = exp(-1 * etzetlam(:,1) * myexp)
            myang(:) = ang_coeff(:) * mya(:) * myb(:) * myfc
            !
            ! Before we multiply by the smoothing functions, we can take
            ! advantage of this calculation and calculate the derivatives
            ! with respect to fc
            !

            ! l is counter for looping over i, j, k
            do l=1, 3
               do m=1, num_etzetlam
                tmp_dang(:,m,l) = mya(j) * myb(m) * mydfc(:,l)
                tmp_dang(:,m,l) = tmp_dang(:,m,l) + &
                    myb(m) * -2 * etzetlam(m,1) * (rij(i,j) + rij(i,k) + rij(j,k)) * mydexp(:,l)
                tmp_dang(:,m,l) = tmp_dang(:,m,l) + &
                    etzetlam(m,2) * etzetlam(m,3) * (1 + etzetlam(m,3) * mycos) ** (etzetlam(m,2) - 1) * mydcos(:,l)
                tmp_dang(:,m,l) = tmp_dang(:,m,l) * ang_coeff(m)
               enddo
            enddo
            !
            ! find the proper angular basis indicies and save the basis
            !
            a = ang_b_ind(my_type, el_key(i))
            b = a + ang_size(my_type, el_key(i)) - 1
            c = atom_basis_index(i)
            e = el_key(i)
            !$OMP CRITICAL
            ang_bas(a:b, c, e) = ang_bas(a:b, c, e) + myang(:)
            !
            ! store the gradient
            !
            ang_grad(:, i, a:b, c, e) = tmp_dang(:,:,1)
            ang_grad(:, j, a:b, c, e) = tmp_dang(:,:,2)
            ang_grad(:, k, a:b, c, e) = tmp_dang(:,:,3)
            !$OMP END CRITICAL

        enddo calc_ang
!$OMP   ENDDO
!$OMP END PARALLEL
        continue

end subroutine calc_bp
!
! Reads the input file
! Calls routines to allocate arrays
! Fills them
subroutine read_input_file(in_path, h5_path)
    ! Read the input file to the bp program and intitialize the variables
    ! Set the requisite variables
    implicit none
    ! I/O Vars
    character(len=:),allocatable, intent(inout) :: in_path, h5_path
    character(len=400) :: line
    integer :: io ! the iostat
    character*200 :: words(6)

    ! Begin
    open(file=in_path, unit=inf, err=1000)
    ! Loop through the file reading it.
    io = 0
    do while (io == 0)
        read(inf,*, iostat=io) line
        if (io .ne. 0) then
            goto 1200
        endif
        call To_lower(line)
        read(line,*) words
        if (words(1) .eq. '#') then
            cycle
        elseif (words(1) .eq. '') then
            cycle
!        elseif (words(1)(1:7) .eq. ('in_file')) then
!            read(words(2), *) h5_file
        elseif (words(1)(1:2) .eq. ('rc')) then
            read(words(2), *, err=1100) Rc
        elseif (words(1)(1:7) .eq. ('num_els')) then
            read(words(2), *, err=1150) num_els
        elseif (words(1)(1:8) .eq. ('max_bond')) then
            read(words(2), *, err=1250) max_bonds
        elseif (words(1)(1:9) .eq. ('max_angle')) then
            read(words(2), *, err=1300) max_angles
        elseif (words(1)(1:10) .eq. ('max_rs_eta')) then
            read(words(2), *, err=1350) max_etas
        elseif (words(1)(1:16) .eq. ('max_zeta_lam_eta')) then
            read(words(2), *, err=1400) max_ezl
        elseif (words(1)(1:) .eq. ('timing_interval')) then
            dotime = .true.
            read(words(2), *, err=1450) timing_interval
        elseif (words(1)(1:7) .eq. ('element')) then
            exit
        else
            goto 1050
        endif
    enddo

    call check_inputs(h5_file)
    call setup_variables()
    call read_elements(ifi)

    close(inf)

    return

   ! Error handling
1000 print *, 'Error opening input file : ', in_path
    stop 'read_input_file 1'

1050 print *, 'Could not read keyword: ', words(1)
    stop 'read_input_file 2'

1100 print*, 'Error reading keyword "rc". Was it a real number?'
    stop 'read_input_file 3'

1150 print*, 'Error reading keyword "num_els". Was it an integer?'
    stop 'read_input_file 4'

1200 print *, 'Unspecified error reading input file. Were there elements in it?'
    stop 'read_input_file 5'

1250 print *, 'Error reading keyword "max_bond". Was it an integer?'
    stop 'read_input_file 6'

1300 print *, 'Error reading keyword "max_angle". Was it an integer?'
    stop 'read_input_file 7'

1350 print *, 'Error reading keyword "max_rs_eta". Was it an integer?'
    stop 'read_input_file 8'

1400 print *, 'Error reading keyword "max_zeta_lam_eta". Was it an integer?'
    stop 'read_input_file 9'

1450 print *, 'Error reading keyword "timing_interval". Was it an integer?'
    stop 'read_input_file 10'

end subroutine read_input_file
!
! Lower-case the string
subroutine To_lower(str)
    character(*), intent(inout) :: str
    integer :: i

    do i = 1, len(str)
        select case(str(i:i))
            case("A":"Z")
            str(i:i) = achar(iachar(str(i:i))+32)
        end select
    end do
end subroutine To_Lower
!
! Validate inputs before going further
subroutine check_inputs(h5_path)
    implicit none
    character(len=:),allocatable, intent(inout) :: h5_path
!    if (.not. allocated(h5_path)) then
!        print *, "Error finding h5 input file, was the keyword in_file not ", &
!                 "given?"
!    endif
    if (num_els .le. 0) then
        print *, "Error, the number of elements must be specified at the top", &
            ' of the file'
        stop 'check_inputs 2'
    endif
    if (max_bonds .le. 0) then
        print *, "Error, the maximum number of bond types must be specified at the top", &
            ' of the file'
        stop 'check_inputs 2'
    endif
    if (max_etas .le. 0) then
        print *, "Error, the maximum number of eta/rs pairs must be specified at the top", &
            ' of the file'
        stop 'check_inputs 2'
    endif
    if (max_angles .le. 0) then
        print *, "Error, the maximum number of angle types must be specified at the top", &
            ' of the file'
        stop 'check_inputs 2'
    endif
    if (max_ezl .le. 0) then
        print *, "Error, the maximum number of zeta/lambda/eta triplets must be specified at the top", &
            ' of the file'
        stop 'check_inputs 2'
    endif
    if (Rc .le. 0) then
        print *, "Error, the Rc cutoff must be specified at the top", &
            ' of the file and be greater than 0'
        stop 'check_inputs 2'
    endif
end subroutine check_inputs
!
! Allocate arrays that we need for reading data from the element section of the
! Input file
!
subroutine setup_variables()
    implicit none
    pi_div_Rc = pi / Rc
    allocate(els(num_els))
    !
    ! Allocate the radial variables
    allocate(eta(max_etas, max_bonds, num_els))
    allocate(rs(max_etas, max_bonds, num_els))
    eta = 0d0
    rs = 0d0
    allocate(rad_size(max_bonds, num_els))
    rad_size = 0
    allocate(radbas_length(num_els))
    radbas_length = 0
    allocate(num_bonds(num_els))
    num_bonds = 0
    allocate(rad_b_ind(max_bonds, num_els))

    !
    ! Allocate the angular variables
    !allocate(zetas(max_zetas, max_angles, num_els))
    !allocate(eprimes(max_eprimes, max_angles, num_els))
    !allocate(lambdas(max_lambdas, max_angles, num_els))
    allocate(ang_size(max_angles, num_els))
    allocate(angbas_length(num_els))
    allocate(ang_types(2, max_angles, num_els))
    allocate(num_angles(num_els))
    allocate(ang_b_ind(max_angles, num_els))
    allocate(etzetlam(max_ezl, 3, max_angles, num_els))
end subroutine setup_variables
!
! Reads the elements section of the input file
subroutine read_elements(ifi)
    implicit none
!
! Local variables
!
    character(len=400) :: line   ! line to read stuff into
    integer :: io                ! the iostat
    character*200 :: words(6)
    logical :: found             ! For exiting Loop
    integer :: el, bt, at, i, re ! Counters

    ! Go back to first element line
    rewind(ifi)
    found = .false.
    do while (io .eq. 0)
        read(ifi, *, err=100, iostat=io) line
        read(line, *) words
        call To_lower(words(1))
        if (words(1)(:7) .eq. 'element') exit
    enddo
    if (io .ne. 0) goto 100

    do el=1, num_el
        !
        ! Begin reading the elements
        read(words(2), *, err=110) els(i)
        read(words(3), *, err=120) num_bonds(el)
        read(words(4), *, err=130) radbas_length(el)
        read(words(5), *, err=140) num_angles(el)
        read(words(6), *, err=150) angbas_length(el)
        !
        ! Read the bonds
        do bt=1, num_bonds(el)
            read(ifi, *, err=160) line
            read(line, *) words
            call to_lower(words(1))
            if (words(1)(1:4) .ne. 'bond') goto 160
            read(words(2), *, err=170) rad_types(bt, el)
            read(ifi, *) line
            read(line, *) words

            re = 0
            do while (words(1) .ne. 'end')
                re = re + 1
                read(words(1), *, err=180) rs(re, bt, el)
                read(words(2), *, err=190) eta(re, bt, el)
                read(ifi, *, iostat=io) line
                if (io .ne. 0) goto 200
                read(line, *) words
                call To_lower(words(1))
            enddo
            rad_size(bt, el) = re

        enddo

        do at=1, num_angle(el)
            read(ifi, *, err=210) line
            read(line, *) words
            call to_lower(words(1))
            if (words(1)(1:4) .ne. 'angle') goto 210
            read(words(2), *, err=220) ang_types(1, at, el)
            read(words(3), *, err=220) ang_types(2, at, el)
            read(ifi, *) line
            read(line, *) words

            re = 0
            do while (words(1) .ne. 'end')
                re = re + 1
                read(words(1), *, err=230) etzetlam(re, 1, at, el)
                read(words(2), *, err=240) etzetlam(re, 2, at, el)
                read(words(3), *, err=250) etzetlam(re, 3, at, el)
                read(ifi, *, iostat=io) line
                if (io .ne. 0) goto 250
                read(line, *) words
                call To_lower(words(1))
            enddo
            ang_size(at, el) = re
        enddo

    enddo



    return
100 print *, 'Error finding elements line on second go. This is more than ', &
             'just a user error'
    stop 'read_elements 1'
110 print *, 'Error on element line ', &
             'First int after element is the atomic number'
    stop 'read_elements 2'
120 print *, 'Error on element line ', &
             'Second int after element is the number of bond types'
    stop 'read_elements 3'
130 print *, 'Error on element line ', &
             'Third int after element is the radial basis length'
    stop 'read_elements 4'
140 print *, 'Error on element line ', &
             'Fourth int after element is the number of angle types'
    stop 'read_elements 5'
150 print *, 'Error on element line ', &
             'Fifth int after element is the angular basis length'
    stop 'read_elements 6'
160 print *, 'Error starting bonds section. Looking for "bond" in line'
    stop 'read_elements 7'
170 print *, 'Error reading atomic number for bond type after "bond" keyword'
    stop 'read_elements 8'
180 print *, 'Error reading rs'
    stop 'read_elements 9'
190 print *, 'Error reading eta'
    stop 'read_elements 10'
200 print *, 'File ended while reading rs and eta during "bond" keyword'
    stop 'read_elements 11'

210 print *, 'Error starting angle section. Looking for "angle" in line'
    stop 'read_elements 12'

215 print *, 'Error reading 2 atomic numbers for angle type after "angle" keyword'
    stop 'read_elements 11'

220 print *, 'File ended while reading zetr, lambda and eta during "angle" keyword'
    stop 'read_elements 11'

230 print *, 'File ended while reading zeta during "angles" keyword'
    stop 'read_elements 11'

240 print *, 'File ended while reading lambda during "angles" keyword'
    stop 'read_elements 11'

250 print *, 'File ended while reading eta during "angles" keyword'
    stop 'read_elements 11'

    end subroutine read_elements

end module bp_symfuncs
