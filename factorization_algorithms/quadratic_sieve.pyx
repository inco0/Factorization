#cython: language_level=3
import logging
from math import ceil, gcd
from random import randint
from gmpy2 cimport *
from gmpy2 import mpz, is_prime, sqrt, exp, log, is_congruent, powmod
from exceptions.exceptions import InvalidInput

logger = logging.getLogger("app")


cpdef list get_quadratic_sieve_factorization(number_to_be_factored: int):
    if number_to_be_factored % 2 == 0 or is_prime(number_to_be_factored):
        raise InvalidInput("Enter an odd composite number that is not a power")

    factors: list[mpz] = []
    cpdef int smooth_numbers_found = 0
    cpdef int required_smooth_numbers = 1
    cpdef int smooth_boundary = 0

    while smooth_numbers_found < required_smooth_numbers:
        logger.info("\n=======================BEGINNING FACTORIZATION=======================\n")
        logger.info("1.PERFORMING INITIALIZATION")
        smooth_boundary, square_roots, factor_base, sieve_length = initialization(number_to_be_factored, smooth_boundary)
        logger.info("2.PERFORMING SIEVING")
        smooth_numbers, smooth_sequence, factorized_smooth_sequence = sieving(smooth_boundary, number_to_be_factored, square_roots,
                                                                sieve_length, len(factor_base))
        smooth_numbers_found = len(smooth_numbers)
        required_smooth_numbers = len(factor_base) + 1
        smooth_boundary *= 4
        logger.debug(f"Found smooth numbers,{smooth_numbers} and smooth sequence factorization {factorized_smooth_sequence}")
    logger.info("3.PERFORMING LINEAR ALGEBRA")
    dependent_rows = linear_algebra(smooth_numbers, factorized_smooth_sequence, factor_base)
    for row_set in dependent_rows:
        factor = factorize(row_set, smooth_numbers, smooth_sequence, number_to_be_factored)
        if factor == 1 or factor == 0 or factor == number_to_be_factored:
            logger.info("Factoring resulted in trivial solution, trying a different combination of rows.")
            pass
        else:
            return [factor, number_to_be_factored//factor]
    logger.info("No set of linearly dependent rows was found.")
    return []


cpdef tuple initialization(number_to_be_factored, smooth_boundary):
    """
    Initializes the variables needed to perform the quadratic sieve factorization

    :param smooth_boundary: If there is a value passed, do not calculate the default value 
    :param number_to_be_factored: The number being factored
    :return: A tuple with smooth boundary B ,a list of roots such that root_i^2≡n(mod p_i),
    the primes up to smooth boundary and the sieve length
    """
    cdef mpz n = mpz(number_to_be_factored)
    if smooth_boundary == 0:
        smooth_boundary = int(ceil(
            sqrt(exp(sqrt(log(n) * log(log(n)))))))  # The value of the smooth boundary B is √(e^(√(ln(n)ln(ln(n))))
    logger.debug(f"The smooth boundary is {smooth_boundary}")
    factor_base = generate_factor_base(number_to_be_factored, smooth_boundary)
    square_roots = get_square_roots(number_to_be_factored, factor_base)
    cpdef long sieve_length = smooth_boundary * 3
    return smooth_boundary, square_roots, factor_base, sieve_length


cpdef tuple sieving(smooth_boundary, number_to_be_factored, square_roots, sieve_length, factor_base_length):
    """
    Sieves the sequence x^2-n for B-smooth values
    
    :param smooth_boundary: The B smooth boundary deciding the upper limit of primes we are sieving with
    :param number_to_be_factored: The number that is being factored
    :param square_roots: A list of tuples with the square roots of n modulo the primes from the prime base
    :param sieve_length: The range of numbers up to which we search B smooth numbers
    :param factor_base_length: The length of the list with the primes up to number_to_be_factored
    :return: A list of B-smooth numbers, and the factors of the smooth numbers sequence
    """
    cpdef mpz initial_sieve_point = mpz(ceil(sqrt(number_to_be_factored)))
    cpdef int smooth_numbers_found = 0
    sieve_sequence = [x * x - number_to_be_factored for x in
                      range(initial_sieve_point, initial_sieve_point + sieve_length)]
    original_sieve_sequence = sieve_sequence.copy()
    sieve_numbers = [x for x in range(initial_sieve_point, initial_sieve_point + sieve_length)]
    prime_factorization = [[] for _ in range(initial_sieve_point, initial_sieve_point + sieve_length)]
    smooth_numbers = []  # Numbers x such that x^2-n is B-smooth
    smooth_sequence = [] # The sequence of numbers x^2-n that are B-smooth
    factorized_smooth_sequence = []  # The numbers x^2-n that are B-smooth factored by the factor base

    # Since we increase the sieve sequence by 1 the result will alternate with an even and odd number
    # So we need to find the first even number and then iterate by 2 through the list to sieve for prime=2
    if sieve_sequence[0] % 2 == 0:
        first_even_index = 0
    else:
        first_even_index = 1
    for i in range(first_even_index, sieve_length, 2):
        while sieve_sequence[i] % 2 == 0:
            sieve_sequence[i] //= 2
            prime_factorization[i].append(2)

    for root, factor in square_roots:
        for i in range((root - initial_sieve_point) % factor, sieve_length, factor):
            while sieve_sequence[i] % factor == 0:
                sieve_sequence[i] //= factor
                prime_factorization[i].append(factor)

    for i in range(sieve_length):
        if sieve_sequence[i] == 1 and prime_factorization[i] != []:
            smooth_numbers.append(sieve_numbers[i])
            factorized_smooth_sequence.append(prime_factorization[i])
            smooth_sequence.append(original_sieve_sequence[i])
            smooth_numbers_found += 1
        if smooth_numbers_found == factor_base_length + 1:
            break

    return smooth_numbers, smooth_sequence, factorized_smooth_sequence


cpdef linear_algebra(smooth_sequence, factorized_smooth_sequence, factor_base):
    """
    Solves the gauss elimination of a matrix with rows being the exponent vectors reduced mod 2 and 
    returns the set of rows that are dependent
    """
    matrix = create_matrix(smooth_sequence, factorized_smooth_sequence, factor_base)
    rows = len(matrix)
    columns = len(matrix[0])
    rows_marked = [False for r in range(rows)]
    for c in range(columns):
        for r in range(rows):
            if matrix[r][c] == 1:
                rows_marked[r] = True
                for k in range(columns):
                    if matrix[r][k] == 1 and k != c:
                        for i in range(rows):
                            matrix[i][k] = (matrix[i][k] + matrix[i][c]) % 2
                break

    dependent_rows = []
    for row, bool_val in enumerate(rows_marked):
        if bool_val == False:
            dependent_rows_set = {row}
            for col in range(columns):
                if matrix[row][col] == 1:
                    for r in range(rows):
                        if matrix[r][col] == 1 and r != row and rows_marked[r] == True:
                            dependent_rows_set.add(r)
            dependent_rows.append(dependent_rows_set)
    return dependent_rows


cpdef create_matrix(smooth_sequence, factorized_smooth_sequence, factor_base):
    """
    Creates and returns a matrix from the prime factorization of the smooth values with
    every row being the exponent vector reduced mod 2
    """
    matrix = []
    for i in range(len(smooth_sequence)):
        exp_vector = []
        for factor in factor_base:
            count = factorized_smooth_sequence[i].count(factor)
            exp_vector.append(count % 2) 
        matrix.append(exp_vector)
        
    return matrix


cpdef legendre(num, p):
    """
    Returns the value of the legendre symbol between n and factor based on the following cases:
    1 if n is a quadratic residue modulo prime and n≢0(mod prime)
    -1 if n is a quadratic non-residue modulo prime
    0 if n ≡ 0(mod prime)
    
    :param num: An integer
    :param p: An odd prime number
    :return: The value of the legendre symbol between n and prime
    """
    cdef prime = mpz(p)
    cdef mpz n = mpz(num % prime)
    cdef mpz t = mpz(1)
    while n != 0:
        while n % 2 == 0:
            n //= 2
            if 3 <= (prime % 8) <= 5:
                t = -t
        n, prime = prime, n
        if is_congruent(n, 3, 4) and is_congruent(prime, 3, 4):
            t = -t
        n = n % prime
    if prime == 1:
        return t
    return 0


cpdef list generate_factor_base(number_to_be_factored, bound):
    """
    Perform the sieve of eratosthenes up to the bound and return all the primes for which the legendre value is 1

    :param number_to_be_factored: The number that is being factored
    :param bound: The number up to which you sieve for primes
    :return: A list with the primes from 2 up to the bound 
    """
    prime_flag = [True] * bound
    primes = []
    for i in range(2, mpz(ceil(sqrt(bound)))):
        if prime_flag[i]:
            for j in range(i * i, bound, i):
                prime_flag[j] = False
    logger.info("========FACTOR BASE FOUND========")
    return [prime for prime in range(2, bound) if
            prime_flag[prime] == True and legendre(number_to_be_factored, prime) == 1]


cpdef list get_square_roots(number_to_be_factored, factor_base):
    """
    Get the two roots of a where a^2 ≡ n(mod prime)

    :param number_to_be_factored: The number that will be factored
    :param factor_base: A list with the primes up to number_to_be_factored
    :return: A list of roots 
    """
    roots = []
    for factor in factor_base[1:]:
        for root in square_root_modulo_prime(number_to_be_factored, factor):
            roots.append(root)
    logger.info("========SQUARE ROOTS FOUND========")
    return roots


cpdef list[tuple] square_root_modulo_prime(number_to_be_factored, prime):
    """
    Given an odd prime and an integer n with (n/p) = 1, this algorithm returns a solution a to a^2 ≡ n (mod prime).
    
    :param number_to_be_factored: The number that is being factored
    :param prime: A prime number, part of the prime base
    :return: The negative and positive solutions of a^2 ≡ n (mod prime) and the prime in a list of tuples
    """
    cdef mpz n = mpz(number_to_be_factored % prime)
    if is_congruent(prime, 3, 8) | is_congruent(prime, 7, 8):
        x = powmod(n, (prime + 1) // 4, prime)
        return [(x, prime), (prime - x % prime, prime)]
    elif is_congruent(prime, 5, 8):
        x = powmod(n, (prime + 3) // 8, prime)
        c = powmod(x, 2, prime)
        if not is_congruent(c, n, prime):
            x = x * powmod(2, (prime - 1) // 4, prime)
        return [(x, prime), (prime - x % prime, prime)]
    else:  # is_congruent(prime, 1, 8)
        while True:
            d = randint(2, prime - 1)
            if legendre(d, prime) == -1:
                break
        s, t = represent_as_power_of_two(prime - 1)
        a = powmod(n, t, prime)
        d = powmod(d, t, prime)
        m = 0
        for i in range(s):
            if is_congruent(int(pow(a * pow(d, m), pow(2, s - i - 1))), -1, prime):
                m += pow(2, i)
        x = pow(mpz(n), mpz((t + 1) // 2)) * pow(mpz(d), mpz(m // 2)) % prime
        return [(x, prime), (prime - x % prime, prime)]


cpdef int factorize(row_set, smooth_numbers, smooth_sequence, number_to_be_factored):
    """
    Calculates the factorization using gcd from the linear dependent rows found from the linear algebra step
    
    :param row_set: The set of rows that are linearly dependent
    :param smooth_numbers: The smooth numbers found from the sieving
    :param smooth_sequence: The value of the sequence of the smooth numbers
    :param number_to_be_factored: The number that is being factored
    :return: A factor of the number that is being factored or 1
    """
    cdef mpz y = mpz(1)
    cdef mpz x = mpz(1)
    for i in row_set:
        y *= smooth_numbers[i]
        x *= smooth_sequence[i]
    x = x % number_to_be_factored
    y = mpz(ceil(sqrt(y)))
    y = y % number_to_be_factored
    
    factor = gcd(x - y, number_to_be_factored)
    return factor


cpdef represent_as_power_of_two(n):
    """
    Represents n as 2^s*t with t odd
    :param n: 
    :return: 
    """
    cpdef int s = 1
    cpdef int t = n // pow(2, s)
    while t % 2 == 0:
        s += 1
        t = n // pow(2, s)
    return s, t
