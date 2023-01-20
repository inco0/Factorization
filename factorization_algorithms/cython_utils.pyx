from gmpy2 cimport *
# import_gmpy2()
#
#
# cdef extern from "gmp.h":
#     void mpz_set_si(mpz_t, long)
#
#
# cpdef get_big_cython_int(n: int):
#     cdef mpz z = GMPy_MPZ_New(NULL)
#     return mpz_set_si(MPZ(z), n)
