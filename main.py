from gmpy2 import is_prime
from exceptions.exceptions import InputIsPrimeNumber
from factorization_algorithms.trial_division import get_trial_division_factorization
from factorization_algorithms.fermat_square_difference import get_fermat_factorization, get_hart_factorization
from factorization_algorithms.pollard_rho import get_pollard_rho_factorization
from factorization_algorithms.quadratic_sieve import get_quadratic_sieve_factorization
from enum import Enum
import logging


class Algorithms(Enum):
    TRIAL_DIVISION = 1
    FERMAT_SQ_DIFF = 2
    HART_SQ_DIFF = 3
    POLLARD_RHO = 4
    QUADRATIC_SIEVE = 5


def trial_division(n: int):
    factors: list = get_trial_division_factorization(n)
    result = 1
    print_factors(n, factors)


def fermat(n: int):
    factors: list = get_fermat_factorization(n)
    result = 1
    print_factors(n, factors)


def hart(n: int):
    factors: list = get_hart_factorization(n, 10000)
    result = 1
    print_factors(n, factors)
    

def pollard_rho(n: int):
    factors = get_pollard_rho_factorization(n)
    print_factors(n, factors)


def quadratic_sieve(n: int):
    factors = get_quadratic_sieve_factorization(n)
    print_factors(n, factors)


def print_factors(number_to_be_factored: int, factors: list):
    print(f"The factors of {number_to_be_factored} are:", end = " ")
    for factor in factors:
        print(factor, end = " ")


def setup_logger():
    logger = logging.getLogger("app")
    logging.basicConfig(format="%(message)s")
    logging.basicConfig()
    logger.setLevel(logging.INFO)


def read_input() -> tuple:
    print(f"Which algorithm would you like to use from the following: \n \
    1. Trial division \n \
    2. Fermat's square difference \n \
    3. Hart's square difference \n \
    4. Pollard Rho's factorization \n \
    5. Quadratic sieve \n ")
    algorithm = int(input("Please enter the number corresponding to the algorithm: "))
    number_to_be_factored = int(input("What number would you like to factor: "))
    if is_prime(number_to_be_factored):
        raise InputIsPrimeNumber()
    return number_to_be_factored, algorithm


def perform_factorization(number_to_be_factored, algorithm):
    if algorithm == Algorithms.TRIAL_DIVISION.value:
        trial_division(number_to_be_factored)
    elif algorithm == Algorithms.FERMAT_SQ_DIFF.value:
        fermat(number_to_be_factored)
    elif algorithm == Algorithms.HART_SQ_DIFF.value:
        hart(number_to_be_factored)
    elif algorithm == Algorithms.POLLARD_RHO.value:
        pollard_rho(number_to_be_factored)
    elif algorithm == Algorithms.QUADRATIC_SIEVE.value:
        quadratic_sieve(number_to_be_factored)


if __name__ == "__main__":
    setup_logger()
    try:
        number_to_be_factored, algorithm = read_input()
        perform_factorization(number_to_be_factored, algorithm)
    except InputIsPrimeNumber:
        print("The number you entered is a prime.")
    except ValueError:
        print("Please enter a number.")
    
