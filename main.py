from factorization_algorithms import trial_division, fermat_square_difference, utils
from factorization_algorithms.pollard_rho import factorize, get_pollard_rho_factorization
from gmpy2 import is_prime
from exceptions.exceptions import InputIsPrimeNumber


if __name__ == "__main__":
    number_to_be_factored = 339126523890
    if is_prime(number_to_be_factored):
        raise InputIsPrimeNumber(f"{number_to_be_factored} is prime.")

    factors = get_pollard_rho_factorization(number_to_be_factored)
    result = 1
    for factor in factors:
        print(factor)
        result *= factor
    assert (result == number_to_be_factored)


