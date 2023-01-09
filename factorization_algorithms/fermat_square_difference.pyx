from math import floor, sqrt
import cython

cpdef factorize(composite_number):
    cdef long long x = floor(sqrt(composite_number))
    cdef int t = 2*x + 1
    cdef int r = x*x - composite_number
    cdef int r_root = int(sqrt(abs(r)))
    while r_root*r_root != r:
        r += t
        t += 2
        r_root = int(sqrt(abs(r)))
    x = (t - 1) // 2
    y = int(sqrt(r))
    print x-y
    print x+y
