#cython: language_level=3
from gmpy2 cimport *
from gmpy2 import mpz, mpz_random, random_state, is_prime, gcd

cpdef list factorize(n: int):
    cdef mpz composite_number = mpz(n)
