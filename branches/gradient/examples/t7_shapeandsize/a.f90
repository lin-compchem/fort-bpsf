program a
implicit none
integer :: mat(2,2) = reshape([1, 2, 3, 4], shape(mat))
integer :: vec(3) = [5, 6, 7]
print *, shape(mat)
print *, rank(mat)
print *, shape(vec) + 1

end program a
