from gmpy2 cimport *
import os
import sys

import_gmpy2()
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))


cdef extern from "gmp.h":
    void mpz_set_si(mpz_t, long)


cdef get_big_cython_int(n: int):
    cdef mpz z = GMPy_MPZ_New(NULL)
    return mpz_set_si(MPZ(z), n)

