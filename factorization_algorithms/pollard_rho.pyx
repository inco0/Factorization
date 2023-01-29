#cython: language_level=3
from gmpy2 cimport *
from utils import get_gmpy_random_state
from math import sqrt
from gmpy2 import mpz, mpz_random, random_state, is_prime, gcd


cpdef list factorize(n: int):
    cdef mpz composite_number = mpz(n)
    cdef mpz b = mpz(1 + mpz_random(get_gmpy_random_state(), n-3))
    cdef mpz s = mpz(mpz_random(get_gmpy_random_state(), n-1))
    cdef mpz c = s
    cdef mpz d = s
    cdef probable_factor = 1

    while probable_factor == 1:
        c = pollard_rho_function(c, b, n)
        d = pollard_rho_function(pollard_rho_function(d, b, n) , b, n)
        probable_factor = gcd(c - d, composite_number)

    if probable_factor < composite_number:
        return [probable_factor]

cpdef mpz pollard_rho_function(x: mpz, b: mpz, n: mpz):
    return mpz((x*x + b % n))
