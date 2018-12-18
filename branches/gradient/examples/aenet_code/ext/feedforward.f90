!-----------------------------------------------------------------------
!      feedforward.f90 - feed-forward artificial neural networks
!-----------------------------------------------------------------------
!+ This file is part of the AENET package.
!+
!+ Copyright (C) 2012-2016 Nongnuch Artrith and Alexander Urban
!+
!+ This program is free software: you can redistribute it and/or modify
!+ it under the terms of the GNU General Public License as published by
!+ the Free Software Foundation, either version 3 of the License, or
!+ (at your option) any later version.
!+
!+ This program is distributed in the hope that it will be useful, but
!+ WITHOUT ANY WARRANTY; without even the implied warranty of
!+ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!+ General Public License for more details.
!+
!+ You should have received a copy of the GNU General Public License
!+ along with this program.  If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------
! 2011-10-18 Alexander Urban (AU)
!-----------------------------------------------------------------------

module feedforward

  implicit none

  public  :: new_Network,            &
             del_Network,            &
             save_Network,           &
             save_Network_ASCII,     &
             load_Network,           &
             load_Network_ASCII,     &
             ff_memsize,             &
             ff_print_info,          &
             ff_get_nweights,        &
             ff_init_weights,        &
             ff_random_init_weights, &
             ff_update_weights,      &
             ff_eval,                &
             ff_deriv,               &
             ff_wderiv,              &
             ff_change_activation,   &
             ff_activate

  private :: print_mat

  !--------------------------------------------------------------------!
  !                            Network type                            !
  !--------------------------------------------------------------------!

  type, public :: Network

     !-----------------------------------------------------------------!
     ! init         : .true., if the instance was initialized          !
     ! memsize      : dynamically allocated words for this instance    !
     ! nlayers      : number of layers in the network, inluding input  !
     !                and output layer                                 !
     ! nnodes(i)    : number of nodes (without bias) in the i-th layer !
     ! nnodes_max   : max. number of nodes in a layer of this network  !
     ! f_a(i)       : activation function type for the i-th layer      !
     ! W(i)         : i-th weight of the graph edges (including bias)  !
     ! iw(i)        : index of the last weight of the i-th layer       !
     !                e.g.: W(iw(2)) --> last weight in the 2nd layer  !
     ! Wsize        : total number of weights; Wsize = size(W)         !
     ! value(i)     : if (evaluated) --> value of the i-th neuron      !
     !                if not         --> undefined                     !
     ! deriv(i)     : if (evaluated) --> derivative of the i-th neuron !
     !                if not         --> undefined                     !
     ! iv(i)        : index of the last node in the i-th layer         !
     ! nvalues      : total number of nodes/neurons in the network     !
     ! evaluated    : .true., if the node values and derivatives have  !
     !                been evaluated                                   !
     ! D(i)         : matrices D_ij of partial derivatives wrt. nodes  !
     ! derivatives  : .true., if the D_ij matrices have been evaluated !
     ! work...      : scratch memory for the computations              !
     !-----------------------------------------------------------------!

     logical                                       :: init = .false.
     integer                                       :: memsize

     integer                                       :: nlayers
     integer,          dimension(:),   allocatable :: nnodes
     integer                                       :: nnodes_max

     integer,          dimension(:),   allocatable :: f_a

     double precision, dimension(:),   allocatable :: W
     integer,          dimension(:),   allocatable :: iw
     integer                                       :: Wsize

     double precision, dimension(:),   allocatable :: value
     double precision, dimension(:),   allocatable :: deriv
     integer,          dimension(:),   allocatable :: iv
     integer                                       :: nvalues
     logical                                       :: evaluated

     double precision, dimension(:),   allocatable :: D
     logical                                       :: derivatives

     double precision, dimension(:),   allocatable :: work
     double precision, dimension(:,:), allocatable :: work2, work3, work4

  end type Network

contains !-------------------------------------------------------------!

  !--------------------------------------------------------------------!
  !                     constructor and destructor                     !
  !--------------------------------------------------------------------!

  function new_Network(arch) result(net)

    implicit none

    integer, dimension(:), intent(in) :: arch
    type(Network)                     :: net

    integer :: ilayer

    net%nlayers = size(arch(:))
    allocate(net%nnodes(net%nlayers),  &
             net%iw(net%nlayers),      &
             net%iv(net%nlayers),      &
             net%f_a(2:net%nlayers)    )

    ! number of nodes per layer:
    do ilayer = 1, net%nlayers
       net%nnodes(ilayer) = arch(ilayer)
    end do
    net%nnodes_max = maxval(net%nnodes(:))

    ! set default activation functions:
    ! (all layers tanh(), linear for the output layer)
    do ilayer = 2, net%nlayers-1
       net%f_a(ilayer) = 1
    end do
    net%f_a(net%nlayers) = 0

    ! indices and size of the weight matrices W_i:
    net%Wsize = 0
    net%iw(1) = 0
    do ilayer = 1, net%nlayers-1
       net%Wsize = net%Wsize + (net%nnodes(ilayer) + 1)*net%nnodes(ilayer+1)
       net%iw(ilayer+1) = net%Wsize
    end do

    ! number of neuron values and index:
    net%nvalues = 0
    net%iv(1)   = 0
    do ilayer = 1, net%nlayers-1
       net%nvalues = net%nvalues + net%nnodes(ilayer) + 1
       net%iv(ilayer+1) = net%nvalues
    end do
    net%nvalues = net%nvalues + net%nnodes(net%nlayers)

    allocate(net%W(net%Wsize),                                 &
             net%D(net%Wsize),                                 &
             net%value(net%nvalues),                           &
             net%deriv(net%nvalues),                           &
             net%work(net%nnodes_max+1),                       &
             net%work2(net%nnodes_max+1,net%nnodes_max+1),     &
             net%work3(net%nnodes_max+1,net%nnodes_max+1),     &
             net%work4(net%nnodes_max+1,(net%nnodes_max+1)**2) )

    net%init        = .true.
    net%evaluated   = .false.
    net%derivatives = .false.

    ! number of words allocated for this object are (at least):
    net%memsize = ff_memsize(net)

  end function new_Network

  !--------------------------------------------------------------------!

  subroutine del_Network(net)

    implicit none

    type(Network), intent(inout) :: net

    if (net%init) then
       deallocate(net%nnodes, net%f_a, net%iw, net%iv,          &
                  net%W, net%D, net%value, net%deriv, net%work, &
                  net%work2, net%work3, net%work4 )
       net%memsize = 0
       net%init    = .false.
    end if

  end subroutine del_Network

  !--------------------------------------------------------------------!
  !                      save / restore networks                       !
  !--------------------------------------------------------------------!

  subroutine save_Network(net, file, unit)

    implicit none

    type(Network),              intent(in) :: net
    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit

    integer            :: u
    integer, parameter :: u_sav = 41

    if (.not. net%init) then
       write(0,*) "Error: network not initialized in `save_Network'."
       stop
    end if

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = u_sav
       open(u, file=trim(file), status='replace', &
            form='unformatted', action='write')
    else
       write(0,*) "Error: neither unit nor file specified in `save_Network'."
       stop
    end if

    write(u) net%nlayers
    write(u) net%nnodes_max
    write(u) net%Wsize
    write(u) net%nvalues

    write(u) net%nnodes(:)
    write(u) net%f_a(:)
    write(u) net%iw(:)
    write(u) net%iv(:)
    write(u) net%W(:)

    if (.not. present(unit)) close(u)

  end subroutine save_Network

  !--------------------------------------------------------------------!

  subroutine save_Network_ASCII(net, file, unit)

    implicit none

    type(Network),              intent(in) :: net
    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit

    integer                     :: u
    integer, parameter          :: u_sav = 41
    integer                     :: i
    character(len=*), parameter :: dfrmt = '(4(1x,ES24.17))'
    character(len=*), parameter :: ifrmt = '(4(1x,I17))'

    if (.not. net%init) then
       write(0,*) "Error: network not initialized in `save_Network'."
       stop
    end if

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = u_sav
       open(u, file=trim(file), status='replace', action='write')
    else
       write(0,*) "Error: neither unit nor file specified in `save_Network_ASCII'."
       stop
    end if

    write(u,*) net%nlayers
    write(u,*) net%nnodes_max
    write(u,*) net%Wsize
    write(u,*) net%nvalues

    write(u,ifrmt) (net%nnodes(i), i=1,net%nlayers)
    write(u,ifrmt) (net%f_a(i), i=2,net%nlayers)
    write(u,ifrmt) (net%iw(i), i=1,net%nlayers)
    write(u,ifrmt) (net%iv(i), i=1,net%nlayers)
    write(u,dfrmt) (net%W(i), i=1,net%Wsize)

    if (.not. present(unit)) close(u)

  end subroutine save_Network_ASCII

  !--------------------------------------------------------------------!

  function load_Network(file, unit) result(net)

    implicit none

    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit
    type(Network)                          :: net

    integer                      :: u
    integer,           parameter :: u_sav = 41

    logical :: fexists

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = u_sav
       inquire(file=trim(file), exist=fexists)
       if (.not. fexists) then
          write(0,*) "Error: file not found in `load_Network': ", &
                     trim(file)
          stop
       end if
       open(u, file=trim(file), status='old', &
            form='unformatted', action='read')
    else
       write(0,*) "Error: neither unit nor file specified in `save_Network'."
       stop
    end if

    read(u) net%nlayers
    read(u) net%nnodes_max
    read(u) net%Wsize
    read(u) net%nvalues

    allocate(net%nnodes(net%nlayers),                          &
             net%f_a(2:net%nlayers),                           &
             net%iw(1:net%nlayers),                            &
             net%iv(1:net%nlayers),                            &
             net%W(net%Wsize),                                 &
             net%D(net%Wsize),                                 &
             net%value(net%nvalues),                           &
             net%deriv(net%nvalues),                           &
             net%work(net%nnodes_max+1),                       &
             net%work2(net%nnodes_max+1,net%nnodes_max+1),     &
             net%work3(net%nnodes_max+1,net%nnodes_max+1),     &
             net%work4(net%nnodes_max+1,(net%nnodes_max+1)**2) )

    read(u) net%nnodes(:)
    read(u) net%f_a(:)
    read(u) net%iw(:)
    read(u) net%iv(:)
    read(u) net%W(:)

    if (.not. present(unit)) close(u)

    net%init        = .true.
    net%evaluated   = .false.
    net%derivatives = .false.

    ! number of words allocated for this object are (at least):
    net%memsize = ff_memsize(net)

  end function load_Network

  !--------------------------------------------------------------------!

  function load_Network_ASCII(file, unit) result(net)

    implicit none

    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit
    type(Network)                          :: net

    integer                      :: u
    integer,           parameter :: u_sav = 41

    integer                     :: i
    character(len=*), parameter :: dfrmt = '(4(1x,ES24.17))'
    character(len=*), parameter :: ifrmt = '(4(1x,I17))'
    logical                     :: fexists

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = u_sav
       inquire(file=trim(file), exist=fexists)
       if (.not. fexists) then
          write(0,*) "Error: file not found in `load_Network': ", &
                     trim(file)
          stop
       end if
       open(u, file=trim(file), status='old', action='read')
    else
       write(0,*) "Error: neither unit nor file specified in `save_Network'."
       stop
    end if

    read(u,*) net%nlayers
    read(u,*) net%nnodes_max
    read(u,*) net%Wsize
    read(u,*) net%nvalues

    allocate(net%nnodes(net%nlayers),                          &
             net%f_a(2:net%nlayers),                           &
             net%iw(1:net%nlayers),                            &
             net%iv(1:net%nlayers),                            &
             net%W(net%Wsize),                                 &
             net%D(net%Wsize),                                 &
             net%value(net%nvalues),                           &
             net%deriv(net%nvalues),                           &
             net%work(net%nnodes_max+1),                       &
             net%work2(net%nnodes_max+1,net%nnodes_max+1),     &
             net%work3(net%nnodes_max+1,net%nnodes_max+1),     &
             net%work4(net%nnodes_max+1,(net%nnodes_max+1)**2) )

    read(u,ifrmt) (net%nnodes(i), i=1,net%nlayers)
    read(u,ifrmt) (net%f_a(i), i=2,net%nlayers)
    read(u,ifrmt) (net%iw(i), i=1,net%nlayers)
    read(u,ifrmt) (net%iv(i), i=1,net%nlayers)
    read(u,dfrmt) (net%W(i), i=1,net%Wsize)

    if (.not. present(unit)) close(u)

    net%init        = .true.
    net%evaluated   = .false.
    net%derivatives = .false.

    ! number of words allocated for this object are (at least):
    net%memsize = ff_memsize(net)

  end function load_Network_ASCII

  !--------------------------------------------------------------------!
  !         size of the allocated memory of a Network instance         !
  !--------------------------------------------------------------------!

  function ff_memsize(net) result(memsize)

    implicit none

    type(Network), intent(in) :: net
    integer                   :: memsize

    if (.not. net%init) then
       memsize = 0
    else
       memsize = (4*net%nlayers - 1)           * 1 &
               + (2*net%Wsize + 2*net%nvalues) * 2 &
               + (net%nnodes_max + 1)          * 2 &
               + 2*(net%nnodes_max+1)**2       * 2 &
               + (net%nnodes_max+1)**3         * 2
    end if

  end function ff_memsize

  !--------------------------------------------------------------------!
  !             print information about a Network instance             !
  !--------------------------------------------------------------------!

  subroutine ff_print_info(net, verbose)

    implicit none

    type(Network),     intent(in) :: net
    logical, optional, intent(in) :: verbose

    character(len=100) :: fname
    integer            :: ilayer
    integer            :: iw1, iw2

    if (.not. net%init) then
       write(*,*) 'The network is not initialized.'
       write(*,*)
       return
    end if

    write(*,'(1x,"Number of layers : ",I3)') net%nlayers
    write(*,*)
    write(*,'(1x,"Number of nodes (without bias) ")')
    write(*,'(1x,"and activation type per layer :")')
    write(*,*)
    write(*,'(5x,I3," : ",I5)') 1, net%nnodes(1)
    do ilayer = 2, net%nlayers
       select case(net%f_a(ilayer))
       case(0)
          fname = 'linear function (linear)'
       case(1)
          fname = 'hyperbolic tangent (tanh)'
       case(2)
          fname = 'logistic function (sigmoid)'
       case(3)
          fname = 'scaled hyperbolic tangent (mtanh)'
       case(4)
          fname = 'scaled hyperbolic tangent + linear twisting (twist)'
       end select
       write(*,'(5x,I3," : ",I5,2x,A)') ilayer, net%nnodes(ilayer), trim(fname)
    end do
    write(*,*)

    write(*,'(1x,"Dynamically allocated words : ",I10," (",F10.2," KB)")') &
         net%memsize, dble(net%memsize*8)/1024.0d0
    write(*,*)

    write(*,'(1x,"Total number of weights (incl. bias) : ",I8)') net%Wsize
    write(*,*)

    verbose1 : if (present(verbose) .and. verbose) then
       write(*,'(1x,"Weight matrices:")')
       write(*,*)

       iw1 = 1
       do ilayer = 1, net%nlayers-1
          iw2 = net%iw(ilayer+1)
          call print_mat(reshape(net%W(iw1:iw2), &
               (/net%nnodes(ilayer+1),net%nnodes(ilayer)+1/) ))
          write(*,*)
       end do
    end if verbose1

  end subroutine ff_print_info

  !--------------------------------------------------------------------!
  !                       weight initialization                        !
  !--------------------------------------------------------------------!

  function ff_get_nweights(net) result(nw)

    implicit none

    type(Network), intent(in) :: net
    integer                   :: nw

    if (.not. net%init) then
       write(0,*) "Error: network not initialized in `get_nweights'."
       stop
    end if

    nw = net%Wsize

  end function ff_get_nweights

  !--------------------------------------------------------------------!

  subroutine ff_init_weights(net, nw, weights)

    implicit none

    type(Network),                   intent(inout) :: net
    integer,                         intent(in)    :: nw
    double precision, dimension(nw), intent(in)    :: weights

    if (.not. net%init) then
       write(0,*) "Error: network not initialized in `init_weights'."
       stop
    end if

    if (nw /= net%Wsize) then
       write(0,*) "Error: wrong number of weights in `init_weights'."
       stop
    end if

    net%W(1:nw) = weights(1:nw)

  end subroutine ff_init_weights

  !--------------------------------------------------------------------!

  subroutine ff_random_init_weights(net)

    implicit none

    type(Network),              intent(inout) :: net

    integer,          dimension(:), allocatable :: iseed
    integer                                     :: i, n, iw
    integer                                     :: clock
    double precision                            :: w, a

    if (.not. net%init) then
       write(0,*) "Error: network not initialized in `random_init_weights'."
       stop
    end if

    call random_seed(size=n)
    allocate(iseed(n))
    call system_clock(count=clock)
    i = 0
    iseed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put=iseed)

    iw = 1
    do i = 1, net%nlayers-1
       n = net%nnodes(i)
       ! select a such that the standard deviation is
       ! 1/sqrt(number of input weights including bias)
       a = 1.0d0/sqrt(dble(n+1))*sqrt(12.0d0)
       do while(iw <= net%iw(i+1))
          call random_number(w)
          net%W(iw) = (w - 0.5d0)*a
          iw = iw + 1
       end do
    end do
    if ((iw - 1) /= net%Wsize) then
       write(0,*) "Error in weight initialization: ", iw-1, net%Wsize
       stop
    end if

    deallocate(iseed)

  end subroutine ff_random_init_weights

  !--------------------------------------------------------------------!

  subroutine ff_update_weights(net, nw, W_upd)

    implicit none

    type(Network),                   intent(inout) :: net
    integer,                         intent(in)    :: nw
    double precision, dimension(nw), intent(in)    :: W_upd

    if (.not. net%init) then
       write(0,*) "Error: network not initialized in `update_weights'."
       stop
    end if

    if (nw /= net%Wsize) then
       write(0,*) "Error: wrong number of weights in `update_weights'."
       stop
    end if

    net%W(:) = net%W(:) + W_upd(:)

  end subroutine ff_update_weights

  !--------------------------------------------------------------------!
  !                                                                    !
  !                 evaluation of the network function                 !
  !                                                                    !
  !--------------------------------------------------------------------!

  subroutine ff_eval(net, nx, x, ny, y)

    implicit none

    type(Network),                   intent(inout) :: net
    integer,                         intent(in)    :: nx
    double precision, dimension(nx), intent(in)    :: x
    integer,                         intent(in)    :: ny
    double precision, dimension(nx), intent(out)   :: y

    integer, dimension(2) :: Wshape
    integer               :: iw1, iw2
    integer               :: iv1, iv2
    integer               :: nnodes1, nnodes2
    integer               :: ilayer

    if (.not. net%init) then
       write(0,*) "Error: network not initialized in `eval'."
       stop
    end if

    if (nx /= net%nnodes(1)) then
       write(0,*) "Error: wrong number of input nodes in `eval': ", nx
       stop
    end if

    if (ny /= net%nnodes(net%nlayers)) then
       write(0,*) "Error: wrong number of output nodes in `eval': ", ny
       stop
    end if

    net%value(1:nx) = x(1:nx)

    iw1 = 1
    iv1 = 1
    nnodes2 = 1 ! just to avoid the gfortran warning
    nnodes1 = net%nnodes(1)
    do ilayer = 1, net%nlayers-1
       iw2 = net%iw(ilayer+1)
       iv2 = net%iv(ilayer+1)

       net%value(iv2) = 1.0d0

       nnodes2 = net%nnodes(ilayer+1)
       Wshape(1:2) = (/ nnodes2, nnodes1+1 /)

       net%work(1:nnodes2) = matmul(         &
            reshape(net%W(iw1:iw2), Wshape), &
            net%value(iv1:iv2)               )

       iv1 = iv2 + 1
       iw1 = iw2 + 1

       call ff_activate(net%f_a(ilayer+1),             &
                        net%work(1:nnodes2),           &
                        net%value(iv1:iv1+nnodes2-1),  &
                        net%deriv(iv1:iv1+nnodes2-1)   )

       nnodes1 = nnodes2
    end do

    if (nnodes1 /= ny) stop 999 ! DEBUG only

    y(1:ny) = net%value(iv1:iv1+nnodes2-1)

    net%evaluated   = .true.
    net%derivatives = .false.

  end subroutine ff_eval

  !--------------------------------------------------------------------!
  !            derivative with respect to the input vector             !
  !--------------------------------------------------------------------!

  subroutine ff_deriv(net, nx, ny, dy)

    implicit none

    type(Network),                      intent(inout) :: net
    integer,                            intent(in)    :: nx, ny
    double precision, dimension(ny,nx), intent(out)   :: dy

    integer, dimension(2) :: Wshape
    integer               :: nnodes0, nnodes1, nnodes2
    integer               :: ilayer, i
    integer               :: iv1, iv2, iw1, iw2

    if (.not. net%init) then
       write(0,*) "Error: network not initialized in `deriv'."
       stop
    end if

    if (.not. net%evaluated) then
       write(0,*) "Error: network not evaluated in `deriv'."
       stop
    end if

    if (nx /= net%nnodes(1)) then
       write(0,*) "Error: wrong number of input nodes in `deriv': ", nx
       stop
    end if

    if (ny /= net%nnodes(net%nlayers)) then
       write(0,*) "Error: wrong number of output nodes in `deriv': ", ny
       stop
    end if

    ! reset all derivatives:
    ! (important to get the zeros for the bias weights)
    net%D(:) = 0.0d0

    ! start with the unity matrix:
    net%work4(:,:) = 0.0d0
    do i = 1, net%nnodes_max+1
       net%work4(i,i) = 1.0d0
    end do

    iw1 = 1
    iv1 = 1
    nnodes2 = 1 ! just to avoid the gfortran warning
    nnodes1 = net%nnodes(1)
    nnodes0 = nnodes1
    layers : do ilayer = 1, net%nlayers-1
       iw2 = net%iw(ilayer+1)
       iv2 = net%iv(ilayer+1)
       nnodes2 = net%nnodes(ilayer+1)

       ! construct diagonal matrix from stored derivatives:
       !
       !       / f'_i1    0      ...   0    \
       !      |   0      f'_i2    0   ...    |
       ! F' = |  ...     ...          ...    |
       !       \  0      ...          f'_iN /
       !
       net%work2(1:nnodes2,1:nnodes2) = 0.0d0
       do i = 1, nnodes2
          net%work2(i,i) = net%deriv(iv2+i)
       end do

       ! smaller matrices --> we don't need the bias weights here:
       Wshape(1:2) = (/ nnodes2, nnodes1 /)

       ! matrix of derivatives:
       !
       !             / d n_i1/d n_(i-1)1  ... d n_i1/d n_(i-1)M \
       ! D = F'*W = |      ...                      ...          |
       !             \ d n_iN/d n_(i-1)1  ... d n_iN/d n_(i-1)M /
       !
       net%work3(1:nnodes2,1:nnodes1) = matmul(      &
            net%work2(1:nnodes2,1:nnodes2),          &
            reshape(net%W(iw1:iw2-nnodes2), Wshape)  )

       ! store the result:
       net%D(iw1:iw2-nnodes2) = reshape(net%work3(1:nnodes2,1:nnodes1), &
                                        (/nnodes2*nnodes1/)             )

       ! combine with the matrix of the previous layer:
       net%work4(1:nnodes2,1:nnodes0) = matmul( &
            net%work3(1:nnodes2,1:nnodes1),     &
            net%work4(1:nnodes1,1:nnodes0)      )

       iv1 = iv2 + 1
       iw1 = iw2 + 1

       nnodes1 = nnodes2
    end do layers

    dy(1:ny,1:nx) = net%work4(1:nnodes2,1:nnodes0)

    net%derivatives = .true.

  end subroutine ff_deriv

  !--------------------------------------------------------------------!
  !     derivative of the output nodes with respect to all weights     !
  !--------------------------------------------------------------------!

  subroutine ff_wderiv(net, nw, ny, dy_dw)

    implicit none

    type(Network),                      intent(inout) :: net
    integer,                            intent(in)    :: nw, ny
    double precision, dimension(ny,nw), intent(out)   :: dy_dw

    integer               :: nnodes1, nnodes2, nnodes3, nnodes12
    integer               :: ilayer, i
    integer               :: iv1, iv2, iw1, iw2
    integer               :: in0, in1, in2

    if (.not. net%init) then
       write(0,*) "Error: network not initialized in `wderiv'."
       stop
    end if

    if (.not. net%evaluated) then
       write(0,*) "Error: network not evaluated in `wderiv'."
       stop
    end if

    if (.not. net%derivatives) then
       write(0,*) "Error: derivatives w.r.t. nodes missing in `wderiv'."
       stop
    end if

    if (nw /= net%Wsize) then
       write(0,*) "Error: wrong number of weights in `wderiv': ", nw
       stop
    end if

    if (ny /= net%nnodes(net%nlayers)) then
       write(0,*) "Error: wrong number of output nodes in `deriv': ", ny
       stop
    end if

    ! start with the unity matrix:
    net%work3(:,:) = 0.0d0
    do i = 1, net%nnodes_max+1
       net%work3(i,i) = 1.0d0
    end do

    iw2 = net%iw(net%nlayers)
    iv2 = net%iv(net%nlayers)
    nnodes2 = net%nnodes(net%nlayers)
    nnodes3 = net%nnodes(net%nlayers)
    layers : do ilayer = net%nlayers-1, 1, -1

       iw1 = net%iw(ilayer) + 1
       iv1 = net%iv(ilayer) + 1

       nnodes1 = net%nnodes(ilayer)

       ! partial derivatives of the nodes in layer 'ilayer+1' with
       ! respect to the weights of the previous layer 'ilayer':
       !
       !  / n'_i1 \
       ! |  n'_i2  |    /                                  \
       ! |  ...    | * | n_(i-1)1   n_(i-1)2  ...  n_(i-1)M |
       !  \ n'_iN /     \                                  /
       !
       !    / d n_i1/d w_i1^(i-1)1    ...   d n_i1/d w_i1^(i-1)M \
       ! = |       ...                             ...            |
       !    \ d n_iN/d w_iN^(i-1)1    ...   d n_iN/d w_iN^(i-1)M /
       !
       ! (note: it's nnodes1+1 because of the bias neurons)
       net%work2(1:nnodes2,1:nnodes1+1) = matmul(                &
            reshape(net%deriv(iv2+1:iv2+nnodes2),(/nnodes2,1/)), &
            reshape(net%value(iv1:iv2),(/1,nnodes1+1/))          )

       ! insert zeros for the vanishing derivatives
       ! d n_ij / d w_kl with (ij) /= (kl)
       nnodes12 = (nnodes1+1)*nnodes2
       net%work4(1:nnodes2,1:nnodes12) = 0.0d0
       in0 = 1
       do in1 = 1, nnodes1+1
       do in2 = 1, nnodes2
          net%work4(in2,in0) = net%work2(in2,in1)
          in0 = in0 + 1
       end do
       end do

       ! combine with derivatives wrt the nodes:
       dy_dw(1:nnodes3,iw1:iw2) = matmul(        &
            net%work3(1:nnodes3,1:nnodes2),      &
            net%work4(1:nnodes2,1:nnodes12)      )

       ! propagate partial derivatives to the next layer:
       net%work3(1:nnodes3,1:nnodes1) = matmul(                      &
            net%work3(1:nnodes3,1:nnodes2),                          &
            reshape(net%D(iw1:iw2-nnodes2), (/ nnodes2, nnodes1 /))  )

       iv2 = iv1 - 1
       iw2 = iw1 - 1

       nnodes2 = nnodes1
    end do layers

  end subroutine ff_wderiv

  !--------------------------------------------------------------------!
  !                                                                    !
  !             change the activation function for a layer             !
  !                                                                    !
  ! Allowed function types are:                                        !
  !                                                                    !
  !    'tanh', 't'    --> hyperbolic tangents                          !
  !    'sigmoid', 's' --> sigmoid function                             !
  !    'linear', 'l'  --> linear function                              !
  !    'mtanh', 'm'   --> modified (scaled) tanh                       !
  !    'twist', 'w'   --> modified (scaled) tanh plus twisting term    !
  !                                                                    !
  !--------------------------------------------------------------------!

  subroutine ff_change_activation(net, ilayer, ftype)

    implicit none

    type(Network),    intent(inout) :: net
    integer,          intent(in)    :: ilayer
    character(len=*), intent(in)    :: ftype

    if (.not. net%init) then
       write(0,*) "Error: network not initialized in `change_activation'."
       stop
    end if

    if ((ilayer < 2) .or. (ilayer > net%nlayers)) then
       write(0,*) "Error: invalid layer number in `change_activation': ", &
                  ilayer
       stop
    end if

    select case(ftype)
    case('linear','l')
       net%f_a(ilayer) = 0
    case('tanh','t')
       net%f_a(ilayer) = 1
    case('sigmoid','s')
       net%f_a(ilayer) = 2
    case('mtanh','m')
       net%f_a(ilayer) = 3
    case('twist','w')
       net%f_a(ilayer) = 4
    case default
       write(0,*) "Error: invalid function type in `change_activation': ", &
                  trim(ftype)
    end select

  end subroutine ff_change_activation

  !--------------------------------------------------------------------!
  !                                                                    !
  !                        activation functions                        !
  !                                                                    !
  ! Function types:                                                    !
  !                                                                    !
  !   0 : linear function f(x) = x                                     !
  !   1 : hyperbolic tangent, y in [-1:1]                              !
  !   2 : sigmoid,            y in [ 0:1]                              !
  !   3 : modified tanh,      y in [-1.7159:1.7159]  f(+/-1) = +/-1    !
  !   4 : tanh & linear twisting term                                  !
  ! [3 & 4 Montavon, Orr, Müller Neural Networks: Tricks of the Trade] !
  !                                                                    !
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  !                 return both: value and derivative                  !
  !--------------------------------------------------------------------!

  elemental subroutine ff_activate(t, x, y, dy)

    implicit none

    integer,          intent(in)  :: t
    double precision, intent(in)  :: x
    double precision, intent(out) :: y
    double precision, intent(out) :: dy

    double precision :: tanhbx
    double precision, parameter :: a = 1.7159d0
    double precision, parameter :: b = 0.666666666666667d0
    double precision, parameter :: c = 0.1d0

    select case(t)
    case(0)
       y  = x
       dy = 1.0
    case(1)
       y  = tanh(x)
       dy = 1.0d0 - y*y
    case(2)
       y  = 1.0d0/(1.0d0 + exp(-x))
       dy = y*(1.0d0 - y)
    case(3)
       tanhbx = tanh(b*x)
       y  = a*tanhbx
       dy = a*(1.0d0 - tanhbx*tanhbx)*b
    case(4)
       tanhbx = tanh(b*x)
       y  = a*tanhbx + c*x
       dy = a*(1.0d0 - tanhbx*tanhbx)*b + c
    case default
       y  = 0.0d0
       dy = 0.0d0
    end select

  end subroutine ff_activate

  !--------------------------------------------------------------------!
  !                                                                    !
  !                some utilities, mostly for debugging                !
  !                                                                    !
  !--------------------------------------------------------------------!

  subroutine print_mat(A)

    implicit none

    double precision, dimension(:,:), intent(in) :: A

    integer                     :: n1, n2
    integer                     :: i1, i2

    double precision            :: val
    double precision, parameter :: EPS = 1.0d-12

    n1 = size(A(:,1))
    n2 = size(A(1,:))

    do i1 = 1, n1
       do i2 = 1, n2
          val = A(i1,i2)
          if (abs(val) < EPS) val = 0.0d0
          write(*, '(1x,ES15.8,2x)', advance='no') val
       end do
       write(*,*)
    end do

  end subroutine print_mat

  !--------------------------------------------------------------------!

!!$  subroutine print_vec(A)
!!$
!!$    implicit none
!!$
!!$    double precision, dimension(:), intent(in) :: A
!!$
!!$    integer                     :: n1
!!$    integer                     :: i1
!!$
!!$    double precision            :: val
!!$    double precision, parameter :: EPS = 1.0d-12
!!$
!!$    n1 = size(A(:))
!!$
!!$    do i1 = 1, n1
!!$       val = A(i1)
!!$       if (abs(val) < EPS) val = 0.0d0
!!$       write(*, '(1x,ES15.8,2x)', advance='no') val
!!$    end do
!!$    write(*,*)
!!$
!!$  end subroutine print_vec

end module feedforward
