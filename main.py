from gmpy2 import is_prime
from exceptions.exceptions import InputIsPrimeNumber
from factorization_algorithms.pollard_rho import get_pollard_rho_factorization
from factorization_algorithms.quadratic_sieve import get_quadratic_sieve_factorization
import logging


def pollard_rho(n: int):
    factors = get_pollard_rho_factorization(n)
    result = 1
    for factor in factors:
        logger.info(factor)
        result *= factor
    assert (result == n)


def quadratic_sieve(n: int):
    get_quadratic_sieve_factorization(n)


def setup_logger():
    logger = logging.getLogger("app")
    logging.basicConfig(format="%(message)s")
    logging.basicConfig()
    logger.setLevel(logging.DEBUG)

def read_input() -> tuple:
    print(f"Which algorithm would you like to use from the following: \n \
    1. Trial division \n \
    2. Fermat's square difference \n \
    3. Dixon's factorization \n \
    4. Pollard Rho's factorization \n \
    5. Quadratic sieve \n ")
    algorithm = int(input("Please enter the number corresponding to the algorithm: "))
    number_to_be_factored = int(input("What number would you like to factor: "))
    if is_prime(number_to_be_factored):
        raise InputIsPrimeNumber()

def perform_factorization(number_to_be_factored, algorithm):
    pass
    
if __name__ == "__main__":
    setup_logger()
    try:
        number_to_be_factored, algorithm = read_input()
        perform_factorization(number_to_be_factored, algorithm)
    except InputIsPrimeNumber:
        print("The number you entered is a prime.")
    except ValueError:
        print("Please enter a number.")
    
