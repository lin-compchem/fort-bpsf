module bp_symfuncs
    implicit none
    integer, parameter :: dp = selected_real_kind(15,300)
    integer, parameter :: num_els = 2
    integer, parameter :: els(2) = [1, 8]
    real(8), parameter :: pi = 4 * atan(1.0d0)
    integer, parameter :: max_bas = 150
    logical, parameter :: calc_grad = .true.
    integer, parameter :: debug = 1
    !TODO: Check for good value for max_bas

    ! General Variables
    real(dp) :: Rc = 8.0
    real(dp) :: pi_div_Rc = 0

    ! Radial variables
    integer, parameter :: max_etas = 24
    integer, parameter :: max_bonds = 2
    real(dp) :: eta(max_etas, max_bonds, num_els) = 0d0
    real(dp) :: rs(max_etas, max_bonds, num_els) = 0d0
    integer :: rad_size(max_bonds, num_els) = 0, radbas_length(num_els) = 0, &
&     rad_types(max_bonds, num_els) = 0, num_bonds(num_els) = 0, &
&     rad_b_ind(max_bonds, num_els)

    ! Angular Variables
    integer, parameter :: max_angles = 3
    integer, parameter :: max_lambdas = 2
    integer, parameter :: max_eprimes = 3
    integer, parameter :: max_zetas = 3
    real(dp) :: zetas(max_zetas, max_angles, num_els) = 0d0
    real(dp) :: eprimes(max_eprimes, max_angles, num_els) = 0d0
    real(dp) :: lambdas(max_lambdas, max_angles, num_els) = 0d0
    real(dp), allocatable :: etzetlam(:,:)
    integer :: num_etzetlam
    integer :: ang_size(max_angles, num_els) =0d0, angbas_length(num_els) = 0d0, &
&              ang_types(2, max_angles, num_els),  num_angles(num_els) = 0,  &
               ang_b_ind(max_angles, num_els)
    
    integer :: num_zetas = max_zetas
    integer :: num_lambdas = max_lambdas
    integer :: num_eprimes = max_eprimes


    ! OPENMP Variables
    ! integer, external :: omp_get_thread_num, omp_get_num_threads
    
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

    ! Timing variables
    logical :: dotime = .true.
    integer :: timing_interval = 100

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
        integer i, j, k, l, m
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
        
        ! Local vars
        ! This is the distance between the ith and jth atom for a given geometry
        real*8, dimension(max_atoms, max_atoms) :: rij
        real*8, dimension(3, max_atoms, max_atoms) :: drijdx
        ! assorted counters
        integer :: i, j, k, l, m, g, x, e, a, b, c, me, my_type, b_ind, &
                   i_btype, j_btype, ut
        ! g : the current geometry
        ! ut : number of elements in upper triangle of matrix with n atoms
        ! i_btype/j_btype : the bond type associated with bond i-j, i.e. the
        !                   me-H will have a bond type of 1, me-O 2, etc
        ! x : counter used to calculate ith and jth index via fancy arithmetic
        ! a, b, c, e : counters used to shorten expressions
        ! mytype: my type of angle, analogous to the btype
        !
        !
        ! Smoothing function array. fc(i,j) = fc(j,i) uses cosine func.
        !  dfcdr == dfc(rij)/dRij
        real*8, dimension(max_atoms, max_atoms) :: fc
        real*8, dimension(3, max_atoms, max_atoms) :: dfcdx
        ! Intermediary variables for calculating angular basis functions
        real*8 :: mycos, myexp, myfc, myang(num_etzetlam), mydcos(3,3)
        real*8 :: mya(num_etzetlam), myb(num_etzetlam), myc(num_etzetlam), ang_coeff(num_etzetlam)
        real*8 :: mydfc(3,3), mydexp(3,3), mydzeta(3,3)
        ! mycos: cos(theta)
        ! myexp: sum of squared dists
        ! myfc:  product of smoothing functions
        ! myang: vector of angles for vectorized computation

        !Basis sets to be reduced in openmp loops
        real*8, dimension(max_bas, max_atoms, num_els) :: tmp_rad_bas, tmp_ang_bas
        real*8, dimension(3, max_atoms, max_bas, max_atoms, num_els) :: tmp_rad_grad, tmp_ang_grad

        ! You can look up an atom index and you will get 1 for H, 2 for O
        integer :: el_key(max_atoms)
        ! Number of atoms for a given element for a geometry
        integer :: g_num_of_els(num_els)
        ! The index in the final basis set for a given atom. I.e. the first oxygen
        ! analyzed will have index of 1, second index of 2, ... then the first hydrogen
        ! will have index of one, etc...
        integer :: atom_basis_index(max_atoms)
        ! Temporary counter
        logical :: el_found
        ! Debugging integers
        integer :: ii, jj, kk

        ! List with all of the triplets of atoms. This is equal to
        ! The upper triangle of a NxN matrix times the number of
        ! atoms
        integer :: cos_inds(3, max_atoms**2 * (max_atoms-1)/2 )
        ! Number of triplets, maximum number of triplets
        integer :: num_cos, max_cos
        ! Holds the radial basis during vectorized computation.
        real :: tmp_rad(max_etas), tmp_drad(3, max_etas, 2), term, tmp_dang(3, num_etzetlam, 3)
        integer :: natm
        ! Counters to keep track of elements location in basis set throughout
        ! multiple geometries.
        integer :: bas_s(num_els), i_end
        integer :: ijkpairs(2,3)
        ! Variables for timing
        real(kind=8) :: t
        integer*8 :: now, clock_rate
        integer*8 :: begin_time, begin_loop

    100 format(I3,'   atomind    ',I3,'   geom')
        max_cos = max_atoms**2 * (max_atoms-1)/2
        pi_div_Rc = pi / Rc
        bas_s(:) = 1

        !!!! Con
        rad_bas(1)%g(:,:,:,:) = 0.d0
        rad_bas(2)%g(:,:,:,:) = 0.d0
        ang_bas(1)%b(:,:) = 0.d0
        ang_bas(2)%b(:,:) = 0.d0
        ang_bas(1)%g(:,:,:,:) = 0.d0
        ang_bas(2)%g(:,:,:,:) = 0.d0

        !!! This should be the same between calculations
        ang_coeff(:) = 2 ** (1 - etzetlam(:,2))

        ! set the diagonals of some arrays to 0
        drijdx(:,:,:) = 0d0

        ! Calculate the vectors and rijs
        call tick(begin_time)
        call tick(begin_loop)
    over_geoms: do g=1, num_geoms
        natm = natoms(g)
        ! Calculate the atom numbers and indicies
        g_num_of_els(:) = 0
        tmp_rad_bas(:, :, :) = 0
        tmp_ang_bas(:, :, :) = 0
        tmp_rad_grad(:,:,:,:,:) = 0
        tmp_ang_grad(:,:,:,:,:) = 0
        ! Calculate the number of upper triangle indices
        ut = natm*(natm-1)/2
        ! Go throgh the atoms in the geometry and find out what element each
        ! atom is. Then tally the number of each element.
 do_natoms: do i =1, natm
                el_found = .false.
                do j=1, num_els
                    if (atmnms(i, g) .eq. els(j)) then
                        g_num_of_els(j) = g_num_of_els(j) + 1
                        el_found = .true.
                        el_key(i) = j
                        atom_basis_index(i) = g_num_of_els(j)
                    endif
                enddo
                if (el_found .eqv. .false.) then
                    print *, 'Error, could not find element in geom'
                    write(*,100) i, g
                    stop 'calc_bas 1'
                endif
           enddo do_natoms
        ! endif

         ! generate the cosine inds. This is an array of all of the unique
         ! angles in the system. This is the key for the threads to use
         ! when the parallel parts start.
        call get_triplets(natm, max_cos, cos_inds, num_cos)

   !$OMP PARALLEL private(x, i, j, k, i_btype, j_btype,  tmp_rad,a,b,c,e) &
   !$OMP& private(my_type, mycos, myexp, myfc, myang) &
   !$OMP& private(term, tmp_drad) &
   !$OMP& reduction(+:tmp_rad_bas) reduction(+:tmp_rad_grad) &
   !$OMP& reduction(+:tmp_ang_bas) reduction(+:tmp_ang_grad)
   !$OMP    DO
           ! Calculate the r-ij vectors and distances and cutoff

calc_bonds: do x=0, ut - 1
            call calc_ij_for_ut(i, j, x, natm)

            ! Distances
            drijdx(:, i, j) = coords(:, j, g) - coords(:, i, g)
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
                if (rad_types(k, el_key(i)) .eq. atmnms(j, g)) then
                    i_btype = k
                    exit
                endif
            enddo
            do k=1, num_bonds(el_key(j))
                if (rad_types(k, el_key(j)) .eq. atmnms(i, g)) then
                    j_btype = k
                    exit
                endif
            enddo

            !Below should be changed if the etas change significantly as it
            !Assumes the etas are the same for all elements to save time
            ! This cannot be reordered because we use intermediary calculation
            ! Results for the variables
            tmp_rad(:) = exp(-1 * eta(:,1,1) * (rij(i,j) - rs(:,1,1))**2)

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
                tmp_rad_bas(a:b, c, e) =  tmp_rad_bas(a:b, c, e) + tmp_rad(:)
                tmp_rad_grad(:, i, a:b, c, e) = tmp_rad_grad(:, i, a:b, c, e) + tmp_drad(:, :, 1)
                tmp_rad_grad(:, j, a:b, c, e) = tmp_rad_grad(:, j, a:b, c, e) + tmp_drad(:, :, 2)
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
                tmp_rad_bas(a:b, c, e) =  tmp_rad_bas(a:b, c, e) + tmp_rad(:)
                tmp_rad_grad(:, i, a:b, c, e) = tmp_rad_grad(:, i, a:b, c, e) - tmp_drad(:, :, 1)
                tmp_rad_grad(:, j, a:b, c, e) = tmp_rad_grad(:, j, a:b, c, e) - tmp_drad(:, :, 2)
            endif
            enddo calc_bonds
    !$OMP ENDDO
    !$OMP DO

calc_ang: do x=1, num_cos
            i = cos_inds(1, x)
            j = cos_inds(2, x)
            k = cos_inds(3, x)

            ! Find which angle type this triplet corresponds to
            my_type = 0
            do m=1, num_angles(el_key(i))
                if (ang_types(1, m, el_key(i)) .eq. atmnms(j, g)) then
                    if (ang_types(2, m, el_key(i)) .eq. atmnms(k, g)) then
                        my_type = m
                    endif
                elseif (ang_types(2, m, el_key(i)) .eq. atmnms(j, g)) then
                    if (ang_types(1, m, el_key(i)) .eq. atmnms(k, g)) then
                        my_type = m
                    endif
                endif
            enddo
            if (my_type .eq. 0) cycle
            ! Calculate the ijk pairs for the derivative calculation.
            ! I.e. in the gaussian filter we need ij and ik for i,
            ! ij and jk for j
!            ! Format: (:,1) is i
!            ! [[j i i
!            !   k k j]]
!            ijkpairs(:,1) = [j,k]
!            ijkpairs(:,2) = [i,k]
!            ijkpairs(;,3) = [j,k]
            !
            ! Calculate the cosine term
            !
            mycos = (rij(i,j)**2 + rij(i,k)**2 - rij(j,k)**2) / &
                    (2.0 * rij(i,j) * rij(i,k))
            ! Derivative of cos with respect to j, k. 2 is j, 3 is k
            term = 1/(rij(i,j) * rij(i,k))
            ! Derivative of cos in the j direction
            mydcos(:,2) = -mycos / rij(i,j) * drijdx(:,j,i)
            mydcos(:,2) = mydcos(:,2) + term * (coords(:,k,g) - coords(:,i,g))
            ! Derivative of cos in the k direction
            mydcos(:,3) = -mycos / rij(i,k) * drijdx(:,k,i)
            mydcos(:,3) = mydcos(:,3) + term * (coords(:,j,g) - coords(:,i,g))
            ! Derivative of cos with respect to i
            mydcos(:,1) = -mycos * (rij(i,k) * drijdx(:,i,j) + rij(i,j)*drijdx(:,i,k))
            mydcos(:,1) = term * (mydcos(:,1) + 2*coords(:,i,g) - coords(:,j,g) - coords(:,k,g))

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
                ! REMOVE THIS LATER!!!!!!!!!!!!!!!
                if (((i .eq. 1) .and. (j .eq. 2)).and.((k.eq.3).and.(g.eq.13))) then
                if (l .eq. 1 .and. m .eq. 1) then
                    ii = i
                    jj = j
                    kk = k
                    print *, 'ds for geom', g
                    print *, 'i',ii,'j',jj,'k',kk
                    print *, 'coords, i then j  then k'
                    print*, coords(:,i,g)
                    print*,coords(:,j,g)
                    print*,coords(:,k,g)
                    print*, 'cos', mycos
                    print*, 'gauss', myb(m)
                    print*, 'fc_prod', myfc
                    print*, 'eta', etzetlam(m,1),'zeta', etzetlam(m,2), 'lam', etzetlam(m,3)
                    print*, 'dzetadtheta'
                    print*, etzetlam(m,2) * etzetlam(m,3) * (1 + etzetlam(m,3) * mycos) ** (etzetlam(m,2) - 1)
                    print*, 'dgauss'
                    print*, myb(m) * -2 * etzetlam(m,1) * (rij(ii,jj) + rij(ii,kk) + rij(jj,kk)) * mydexp(:,l)
                    print*, 'dfc'
                    print*,  mya(j) * myb(m) * mydfc(:,l)
                    print*, 'dcosd Xi'
                    print*, mydcos(:,1)
                    print*, 'dcosd Xj'
                    print*, mydcos(:,2)
                    print*, 'dcosd Xk'
                    print*, mydcos(:,3)
                    print*, 'rij', 'rik', 'rjk'
                    print*, rij(ii,jj),rij(ii,kk),rij(jj,kk)
                    print*, 'drijdi'
                    print*, drijdx(:,ii,jj)
                    print*, 'drijdj'
                    print*, drijdx(:,jj,ii)
                    print*, 'drikdi'
                    print*, drijdx(:,ii,kk)
                    print*, 'drikdk'
                    print*, drijdx(:,kk,ii)
                    print*, 'drjkdj'
                    print*, drijdx(:,jj,kk)
                    print*, 'drjkdk'
                    print*, drijdx(:,kk,jj)
                    print*, 'xj-xi or vij'
                    print*, coords(:,j,g) - coords(:,i,g)
                    print*,'xk-xi or vik'
                    print*, coords(:,k,g) - coords(:,i,g)
                    print*,'xk-xj or vjk'
                    print*, coords(:,k,g) - coords(:,j,g)
                endif
                endif
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
            tmp_ang_bas(a:b, c, e) = tmp_ang_bas(a:b, c, e) + myang(:)
            !
            ! store the gradient
            !
            tmp_ang_grad(:, i, a:b, c, e) = tmp_dang(:,:,1)
            tmp_ang_grad(:, j, a:b, c, e) = tmp_dang(:,:,2)
            tmp_ang_grad(:, k, a:b, c, e) = tmp_dang(:,:,3)

        enddo calc_ang
!$OMP   ENDDO
!$OMP END PARALLEL
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
1000 print *, "Error, could not match element"
    goto 9000
9000 print *, i, j, g
    stop 'final_exit'
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
end module bp_symfuncs
