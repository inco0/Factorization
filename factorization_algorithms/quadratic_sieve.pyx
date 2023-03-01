#cython: language_level=3
from gmpy2 cimport *
from gmpy2 import mpz, is_prime, sqrt, ceil, exp, log, is_congruent, powmod, mpz_random
from exceptions.exceptions import InvalidInput
from math import pow
from random import random


cpdef list get_quadratic_sieve_factorization(number_to_be_factored: int):
    factors: list[mpz] = []
    if number_to_be_factored % 2 == 0 or is_prime(number_to_be_factored):
        raise InvalidInput("Enter an odd composite number that is not a power")
    smooth_boundary, roots = initialization(number_to_be_factored)
    # sieving()
    # linear_algebra()
    # factorize()


cpdef tuple initialization(number_to_be_factored):
    """
    Initializes the variables needed to perform the quadratic sieve factorization

    :param number_to_be_factored: The number being factored
    :return: A tuple with B and a list of roots such that root_i^2≡n(mod p_i)
    """
    cdef mpz n = mpz(number_to_be_factored)
    cdef mpz smooth_boundary = mpz(
        ceil(sqrt(exp(sqrt(log(n) * log(log(n)))))))  # The value of B is √(e^(√(ln(n)ln(ln(n))))
    quadratic_roots = get_quadratic_roots(number_to_be_factored, smooth_boundary)

cpdef legendre(num, prime):
    """
    Returns the value of the legendre symbol between n and prime based on the following cases:
    1 if n is a quadratic residue modulo prime and n≢0(mod prime)
    -1 if n is a quadratic non-residue modulo prime
    0 if n ≡ 0(mod prime)
    
    :param num: An integer
    :param prime: An odd prime number
    :return: The value of the legendre symbol between n and prime
    """
    cdef mpz n = mpz(num % prime)
    cdef mpz t = 1
    while n != 0:
        while n % 2 == 0:
            n /= 2
            if 3 <= (prime % 8) <= 5:
                t = -t
        n, prime = prime, n
        if is_congruent(n, 3, 4) and is_congruent(prime, 3, 4):
            t = -t
        n = n % prime
    if prime == 1:
        return t
    return 0

cpdef list sieve_of_eratosthenes(n):
    """
    Perform a sieving up to n and return a list of primes

    :param n: The number up to which you sieve for primes
    :return: A list with the primes from 3 up to n 
    """
    prime_flag = [True] * n
    primes = []
    for i in range(2, ceil(sqrt(n))):
        if prime_flag[i]:
            for j in range(i * i, n, i):
                prime_flag[j] = False
    return [prime for prime in range(n) if prime_flag[prime] == True and prime > 2]

cpdef list get_quadratic_roots(number_to_be_factored, smooth_bound):
    """
    Get the two roots of a where a^2 ≡ n(mod prime)

    :param number_to_be_factored: The number that will be factored
    :param smooth_bound: The B smooth bound
    :return: A list of roots 
    """
    primes = sieve_of_eratosthenes(number_to_be_factored)
    roots = []
    for prime in primes:
        if legendre(number_to_be_factored, prime) == 1:
            pass
    return roots


cpdef square_roots(num, prime):
    """
    Given an odd prime and an integer a with (n/p) = 1, this algorithm returns a solution a to a^2 ≡ n (mod prime).
    Algorithm described in "Prime Numbers, a computational Perspective - Richard Crandall, Carl Pomerance", 2.3.8

    :param num: Finds a such that a^2 ≡ n(mod p)
    :param prime: A prime number part of the prime base
    :return: The negative and positive solutions of a^2 ≡ n (mod prime)
    """
    cdef mpz n = mpz(num % prime)
    if is_congruent(prime, 3, 8) | is_congruent(prime, 7, 8):
        x = powmod(n, (prime + 1) // 4, prime)
        return x, prime - x % prime
    elif is_congruent(prime, 5, 8):
        x = powmod(n, (prime + 3) // 8, prime)
        c = (x*x) % prime
        if not is_congruent(c, n, prime):
            x = (x * powmod(2, (prime - 1) // 4), prime)
        return x, prime - x % prime
    else : # is_congruent(prime, 1, 8)
        while True:
            d = random.randint(2, prime - 1)
            if legendre(d, prime) == -1:
                break
        s, t = represent_as_power_of_two(prime - 1)
        a = powmod(n, t, prime)
        d = powmod(d, t, prime)
        m = 0
        for i in range(s):
            if is_congruent(pow(a*pow(d, m), pow(2, s-i-1)), -1, prime):
                m += pow(2, i)
        x = (pow(n, (t + 1) // 2) * pow(d, m // 2)) % prime
        return x, prime - x % prime
