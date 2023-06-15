#cython: language_level=3
from gmpy2 cimport *
from gmpy2 import random_state
from math import sqrt
from gmpy2 import mpz, mpz_random, random_state, is_prime, gcd
from exceptions.exceptions import FactorizationFailed, UnfinishedFactorization


cpdef int factorize(n: int):
    cdef mpz composite_number = mpz(n)
    cdef mpz b = mpz(1 + mpz_random(random_state(hash(random_state())), n-3))
    cdef mpz s = mpz(mpz_random(random_state(hash(random_state())), n-1))
    cdef mpz c = s
    cdef mpz d = s
    cdef probable_factor = 1

    while probable_factor == 1:
        c = mpz(pollard_rho_function(c, b, composite_number))
        d = mpz(pollard_rho_function(pollard_rho_function(d, b, composite_number) , b, composite_number))
        probable_factor = gcd(c - d, composite_number)

    if probable_factor < composite_number:
        return probable_factor
    else:
        raise UnfinishedFactorization()


cpdef mpz pollard_rho_function(x: mpz, b: mpz, composite_number: mpz):
    return mpz(x*x + b % composite_number)


cpdef list get_pollard_rho_factorization(number_to_be_factored: int):
    cdef mpz i = mpz(0)
    factors: list[mpz] = []
    while i <= 10:
        try:
            factor = factorize(number_to_be_factored)
            factors.append(factor)
            # while is_prime(factor):
            #     factors.append(factor)
            number_to_be_factored //= factor
            if is_prime(number_to_be_factored):
                factors.append(number_to_be_factored)
                return factors
        except UnfinishedFactorization:
            factors: list = factorize(number_to_be_factored)
        i += 1
    raise FactorizationFailed(f"Could not find enough proper factors of {number_to_be_factored}")
