program loop_fuse
  implicit none
  integer i, j, k, n, x, z, ut
  n = 4
  ut = n*(n-1)/2
  print *, "ut iterations", ut
  print *, "total iterations", n*ut
  do z=0, n * ut - 1
     i = z/ut
     x = (z - i * ut) 
     j = x /n
     k = mod(x, n)
     if (k .le. j) then
        j = n - j - 2
        k = n - k - 1
     endif
     i = i + 1
     j = j + 1
     k = k + 1
     print *,"ijkx", i,  j, k, x
  enddo
end program
