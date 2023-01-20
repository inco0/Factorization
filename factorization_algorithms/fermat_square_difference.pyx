#cython: language_level=3
from math import floor, sqrt
from gmpy2 cimport *
import_gmpy2()


cdef extern from "gmp.h":
    void mpz_init(mpz_t)
    void mpz_init_set_si(mpz_t, long)
    void mpz_add(mpz_t, mpz_t, mpz_t)
    void mpz_sub(mpz_t, mpz_t, mpz_t)
    void mpz_mul(mpz_t, mpz_t, mpz_t)
    void mpz_mul_si(mpz_t, mpz_t, long)


cpdef factorize(n: int):
    cdef mpz composite_number = mpz_init_set_si(MPZ(GMPy_MPZ_New(NULL)),
                                                n)
    cdef mpz x = GMPy_MPZ_New(NULL)
    cdef mpz r = GMPy_MPZ_New(NULL)
    cdef mpz t = GMPy_MPZ_New(NULL)
    cdef mpz r_root = GMPy_MPZ_New(NULL)
    print(x)
    # mpz_init_set_si(MPZ(x), floor(sqrt(composite_number)))
    # t = 2*x + 1
    # mpz_init_set_si(mpz_sub(mpz_mul(MPZ(x),x,x), composite_number))
    # r_root = sqrt(abs(r))
    #
    # while r_root*r_root != r:
    #     r += t
    #     t += 2
    #     r_root = int(sqrt(abs(r)))
    # x = (t - 1) // 2
    # y = int(sqrt(r))
    # print(int(x-y))
    # print(x+y)
