from math import floor, sqrt
from gmpy2 cimport *
from cython_utils import get_big_cython_int
import cython

cpdef factorize(composite_number):
    cdef mpz x = get_big_cython_int(floor(sqrt(composite_number)))
    cdef int t = 2*x + 1
    cdef mpz r = get_big_cython_int(x*x - composite_number)
        # x*x - composite_number
    cdef int r_root = int(sqrt(abs(r)))
    while r_root*r_root != r:
        r += t
        t += 2
        r_root = int(sqrt(abs(r)))
    x = (t - 1) // 2
    y = int(sqrt(r))
    print x-y
    print x+y
