from gmpy2 import is_prime
from exceptions.exceptions import InputIsPrimeNumber
from factorization_algorithms.pollard_rho import get_pollard_rho_factorization
from factorization_algorithms.quadratic_sieve import get_quadratic_sieve_factorization


def pollard_rho(n: int):
    factors = get_pollard_rho_factorization(n)
    result = 1
    for factor in factors:
        print(factor)
        result *= factor
    assert (result == n)


def quadratic_sieve(n: int):
    get_quadratic_sieve_factorization(n)


if __name__ == "__main__":
    number_to_be_factored = 800000987987698796780018000009879876987967800180000098798769879678001
    if is_prime(number_to_be_factored):
        raise InputIsPrimeNumber(f"{number_to_be_factored} is prime.")
    # pollard_rho(number_to_be_factored)
    print(f"Initializing quadratic sieve factoring of {number_to_be_factored}")
    quadratic_sieve(number_to_be_factored)
