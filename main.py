from factorization_algorithms import trial_division, fermat_square_difference, utils
from gmpy2 import is_prime


if __name__ == "__main__":
    number_to_be_factorized = 105
    if is_prime(number_to_be_factorized):
        print("Number is most likely prime")
        exit()
    factors: list = trial_division.factorize(number_to_be_factorized)
    result = 1
    for factor in factors:
        print(factor)
        result *= factor
    assert (result == number_to_be_factorized)
