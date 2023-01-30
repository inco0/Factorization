from factorization_algorithms import trial_division, fermat_square_difference, utils
from factorization_algorithms.pollard_rho import factorize, get_polard_rho_factorization
from gmpy2 import is_prime


if __name__ == "__main__":
    number_to_be_factored = 105
    if is_prime(number_to_be_factored):
        print("Number is most likely prime")
        exit()

    get_pollard_rdho_factorization(number_to_be_factored)
    result = 1
    for factor in factors:
        print(factor)
        result *= factor
    assert (result == number_to_be_factored)


