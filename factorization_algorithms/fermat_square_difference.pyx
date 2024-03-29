#cython: language_level=3
from gmpy2 cimport *
from gmpy2 import mpz, is_square, ceil, floor, sqrt, gcd, is_prime


cpdef list get_fermat_factorization(n: int):
    cdef mpz composite_number = mpz(n)
    cdef mpz first_factor = mpz(floor(sqrt(composite_number)))
    cdef mpz r = mpz(first_factor*first_factor - composite_number)
    cdef mpz t = mpz(2*first_factor + 1)
    cdef mpz r_root = mpz(sqrt(abs(r)))
    while r_root*r_root != r:
        print(r)
        r += t
        t += 2
        r_root = mpz(int(sqrt(abs(r))))
    first_factor = (t - 1) // 2
    cdef mpz second_factor = mpz(int(sqrt(r)))
    return [first_factor - second_factor, first_factor + second_factor]

cpdef list get_hart_factorization(n: int, limit: int):
    cdef mpz composite_number = mpz(n)
    cdef mpz s = mpz(0)
    cdef mpz m = mpz(0)
    for i in range(1, limit):
        s = mpz(ceil(sqrt(composite_number*i)))
        m = mpz(s * s % composite_number)
        if is_square(m):
            break
    cdef mpz t = mpz(sqrt(m))
    cdef factor = mpz(gcd(s-t, composite_number))
    return [factor, composite_number // factor]
