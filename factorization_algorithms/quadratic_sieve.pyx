#cython: language_level=3
from gmpy2 cimport *
from gmpy2 import mpz, mpz_random, random_state, is_prime, gcd, sqrt, ceil, exp, log
from exceptions.exceptions import InvalidInput
import numpy as np

cpdef list prime_sieve(smooth_bound)
    """
    """

cpdef tuple initialization(number_to_be_factored):
    """
    Initializes the things needed to do the factorization

    :param number_to_be_factored : The number being factored
    :return : A tuple with B and a list of roots such that root_i^2≡n(mod p_i)
    """
    cdef mpz n = mpz(number_to_be_factored)
    cdef mpz smooth_boundary = mpz(ceil(sqrt(exp(sqrt(log(n)*log(log(n))))))) # The value of B is √(e^(√(ln(n)ln(ln(n))))
    primes = prime_sieve(smooth_boundary)

cpdef list get_quadratic_sieve_factorization(number_to_be_factored: int):
    factors: list[mpz] = []
    if number_to_be_factored % 2 == 0 or is_prime(number_to_be_factored):
        raise InvalidInput("Enter an odd composite number that is not a power")
    b, roots = initialization(number_to_be_factored)
    sieving()
    linear_algebra()
    factorize()
