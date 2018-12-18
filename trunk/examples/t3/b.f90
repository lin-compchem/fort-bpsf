program loop_fuse
  implicit none
  integer i, j, n, x
  n = 4
  do x=0, n*(n-1)/2 - 1
     i = x/n
     j = mod(x, n)
     if (j .le. i) then
        i = n - i - 2
        j = n - j - 1
     endif
     i = i + 1
     j = j + 1
     print *, i, j
  enddo
end program
