from math import sqrt
from random import randint

import factorization.utils as utils
import factorization.primality_tests as primality


def recursive_step(divisor, factorization, primality_tests):
    # Get the factorization of the mcd recursively.
    mcd_factorization = factorize(divisor, primality_tests=primality_tests)
    # append it to the factorization of n.
    factorization = utils.add_divisor_factorization(factorization, mcd_factorization)

    # we do the same with the current reduced value of the factorization.
    reduced_value_factorization = factorize(factorization.reduced_value, primality_tests=primality_tests)
    factorization = utils.add_divisor_factorization(factorization, reduced_value_factorization)

    # We now break the loop and return the factorized value.
    return factorization


def factorize(n, primality_tests=None,
              threshold=None, a=None,
              max_iters=10000):
    """
    Applies the p -1 pollard factorization method
    :param n: Integer to be factorized.
    :param primality_tests: Primality tests to be used.
    :param threshold: initial threshold of the algorithm.
    :param a: initial number of algorithm.
    :param max_iters: maximum number of times the algorithm should be tried.
    :return: A Factorization object with the factorized integer.
    """
    # Initialize factorization object.
    factorization = utils.Factorization(n)
    n = abs(n)  # Set n as a positive number to avoid problems.

    # Initialize primality tests.
    if primality_tests is None:
        primality_tests = [primality.miller_rabin_primality_test]

    # Check if the Factorization object is already decomposed (if its reduced value is prime).
    # If the reduced value is indeed prime we are done,
    # add the reduced value to the factorization and return the factorization.
    if primality.is_decomposed(factorization, primality_tests):
        factorization.add_factor(factorization.reduced_value)
        return factorization

    # If a threshold is not set we set a random threshold ourselves.
    if threshold is None:
        threshold = randint(1, 2 * int(sqrt(n)))

    # If a value a is not set we set a random a value ourselves.
    if a is None:
        # if n is less or equal 3 then it is already factorized.
        if n <= 3:
            factorization.add_factor(n)
            return factorization
        a = randint(2, n - 2)

    # Try the algorithm at most max_iters times.
    for i in range(max_iters):
        # check that a is not itself a divisor of n.
        mcd = utils.euclidean_algorithm_m_c_d(a, n)
        if mcd not in [1, n]:
            # We have found a divisor of n by sheer luck.
            return recursive_step(mcd, factorization, primality_tests)

        # If we are not lucky then we get k = threshold! mod n-1 (We are going to make the k power of a so this
        # is equivalent and will reduce computation time.
        k = utils.get_factorial_modulo(threshold, n - 1)
        # Get a^k -1 modulo n.
        ak = utils.get_power_mod_n(a, k, n) - 1

        # check if ak-1 has common divisors with n.
        mcd = utils.euclidean_algorithm_m_c_d(ak-1, n)
        if mcd not in [1, n]:
            # We have found a divisor of n by sheer luck.
            return recursive_step(mcd, factorization, primality_tests)

        # if we have found nothing we prepare for next iteration.
        # if n is less or equal 3 then it is already factorized.
        if n <= 3:
            factorization.add_factor(n)
            return factorization
        threshold = randint(1, int(sqrt(n)) + 1)  # set a random threshold.
        a = randint(2, n - 2)  # set a random a value.

    # If we get no divisor then we return the factorization as it is with a warning that a
    # different rho function may be able to complete the factorization.
    print('Warning unable to decompose value {} using -1 pollard factorization.'.format(factorization.reduced_value))
    return factorization
