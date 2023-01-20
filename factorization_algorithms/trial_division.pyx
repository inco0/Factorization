#cython: language_level=3
from math import sqrt
from gmpy2 import mpz
import cython

cpdef factorize(composite_number):
    # import pydevd_pycharm
    # pydevd_pycharm.settrace('localhost', port=40050, stdoutToServer=True, stderrToServer=True)
    current_number = composite_number
    cdef int prime = 2
    divisors: [int] = []

    while prime <= sqrt(composite_number):
        if current_number % prime == 0:
            divisors.append(prime)
            current_number //= prime
        else:
            prime += 1
    if current_number == composite_number:
        divisors.append(current_number)
    elif current_number > 1:
        divisors.append(current_number)

    print(f'The divisors of {composite_number} are:', end=' ')
    for divisor in divisors:
        print(f'Divisor', end= ' ')