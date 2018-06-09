from math import sqrt
from random import randint
from typing import List, Callable, Optional

from factorization.utils import Factorization, euclidean_algorithm_m_c_d, add_divisor_factorization
from factorization.primality_tests import miller_rabin_primality_test, is_decomposed, get_power_mod_n


def recursive_step(divisor, factorization: Factorization,
                   primality_tests: List[Callable[[int], bool]]) -> Factorization:
    # Get the factorization of the mcd recursively.
    mcd_factorization = factorize(divisor, primality_tests=primality_tests)
    # append it to the factorization of n.
    factorization = add_divisor_factorization(factorization, mcd_factorization)

    # we do the same with the current reduced value of the factorization.
    reduced_value_factorization = factorize(factorization.reduced_value, primality_tests=primality_tests)
    factorization = add_divisor_factorization(factorization, reduced_value_factorization)

    # We now break the loop and return the factorized value.
    return factorization


def get_factorial_modulo_n(a: int, n: int) -> int:
    prod = 1
    for i in range(1, a + 1):
        prod *= i
        prod = prod % n
    return prod


def factorize(n, primality_tests: Optional[List[Callable[[int], bool]]] = None,
              threshold: Optional[int]=None, a: Optional[int]=None,
              max_iters: int=100) -> Factorization:
    # Initialize factorization object.
    factorization = Factorization(n)
    n = abs(n)  # Set n as a positive number to avoid problems.

    # Initialize primality tests.
    if primality_tests is None:
        primality_tests = [miller_rabin_primality_test]

    # Check if the Factorization object is already decomposed (if its reduced value is prime).
    # If the reduced value is indeed prime we are done,
    # add the reduced value to the factorization and return the factorization.
    if is_decomposed(factorization, primality_tests):
        factorization.add_factor(factorization.reduced_value)
        return factorization

    # If a threshold is not set we set a random threshold ourselves.
    if threshold is None:
        threshold = randint(1, 2*int(sqrt(n)))

    # If a value a is not set we set a random a value ourselves.
    if a is None:
        a = randint(2, n-2)

    # Try the algorithm at most max_iters times.
    for i in range(max_iters):
        # check that a is not itself a divisor of n.
        mcd = euclidean_algorithm_m_c_d(a, n)
        if mcd not in [1, n]:
            # We have found a divisor of n by sheer luck.
            return recursive_step(mcd, factorization, primality_tests)

        # If we are not lucky then we get k = threshold! mod n-1 (We are going to make the k power of a so this
        # is equivalent and will reduce computation time.
        k = get_factorial_modulo_n(threshold, n-1)
        # Get a^k -1 modulo n.
        ak = get_power_mod_n(a, k, n) - 1

        # check if ak-1 has common divisors with n.
        if ak not in [1, n]:
            # We have found a divisor of n by sheer luck.
            return recursive_step(ak, factorization, primality_tests)

        # if we have found nothing we prepare for next iteration.
        threshold = randint(1, int(sqrt(n)) + 1)  # set a random threshold.
        a = randint(2, n - 2)  # set a random a value.

    # If we get no divisor then we return the factorization as it is with a warning that a
    # different rho function may be able to complete the factorization.
    print('Warning unable to decompose value {} using'
          'rho -1 pollard factorization.'.format(factorization.reduced_value))
    return factorization
