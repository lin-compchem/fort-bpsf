program loop_fuse
  implicit none
  integer i, j, n, x
  n = 4
  do x=1, n*n
     i = x/n
     j = mod(x, n)
     print *, i, j
  enddo
end program
