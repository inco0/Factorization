#cython: language_level=3
from gmpy2 import mpz, sqrt
from gmpy2 cimport *


cpdef get_trial_division_factorization(composite_number: int):
    cdef mpz current_number = mpz(composite_number)
    cdef mpz prime = mpz(2)
    divisors: [mpz] = []

    while prime <= mpz(sqrt(composite_number)):
        if current_number % prime == 0:
            divisors.append(prime)
            current_number //= prime
        else:
            prime += 1
    if current_number == composite_number:
        divisors.append(current_number)
    elif current_number > 1:
        divisors.append(current_number)

    return divisors