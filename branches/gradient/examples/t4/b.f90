program loop_fuse
  implicit none
  integer i, j, n, x, t
  n = 4
  t = 16

  do i=1, t
     j = i / 4
     print *, i,j
  enddo

end program
