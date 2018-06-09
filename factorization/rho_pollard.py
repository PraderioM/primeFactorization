from typing import List, Callable, Optional
from random import randint

from factorization.utils import Factorization, euclidean_algorithm_m_c_d, add_divisor_factorization
from factorization.primality_tests import miller_rabin_primality_test, is_decomposed


def factorize(n, primality_tests: Optional[List[Callable[[int], bool]]] = None,
              rho_function: Optional[Callable[[int], int]] = None,
              factorization: Optional[Factorization] = None):
    # Initialize primality tests.
    if primality_tests is None:
        primality_tests = [miller_rabin_primality_test]

    # Initialize rho_function.
    if rho_function is None:
        # Random integer for rho function
        aux = randint(0, n - 1)

        def rho_function(x):
            return (x ** 2 + aux) % n

    # Initialize factorization object. We need to do this to call the function recursively.
    if factorization is None:
        factorization = Factorization(n)

    # Check if the Factorization object is already decomposed (if its reduced value is prime).
    # If the reduced value is indeed prime we are done,
    # add the reduced value to the factorization and return the factorization.
    if is_decomposed(factorization, primality_tests):
        factorization.add_factor(factorization.reduced_value)
        return factorization

    # If not we run a step of the rho Pollard factorization algorithm.
    # We initialize the set of successive images of the rho function with a random integer between 0 and n.
    last_rho_function_image = randint(0, n - 1)
    compare_value = last_rho_function_image
    n_bits = 1  # number of iterations to be performed before next compare value change.
    # The loop cycle can at most be as long as n.
    # If we have to do loops that long we will know wer are stuck in a cycle and will stop.
    # We could find earlier that we are stuck in a cycle by checking if a certain value repeats
    # itself but it would be more time consuming.
    while n_bits < n:
        # We iterate over the n_bits +1 bits digits.
        for i in range(n_bits):
            # get next image of rho function
            last_rho_function_image = rho_function(last_rho_function_image)

            # get maximum common divisor between n and last_rho_function_image-compare_value.
            mcd = euclidean_algorithm_m_c_d(n, last_rho_function_image-compare_value)
            # if the mcd is not trivial we continue have found a divisor and we continue the algorithm recursively.
            if mcd not in [1, n]:
                # Get the factorization of the mcd recursively.
                mcd_factorization = factorize(mcd)
                # append it to the factorization of n.
                factorization = add_divisor_factorization(factorization, mcd_factorization)

                # we do the same with the current reduced value of the factorization.
                reduced_value_factorization = factorize(mcd)
                factorization = add_divisor_factorization(factorization, reduced_value_factorization)

                # We now break the loop and return the factorized value.
                return factorization

            # Else we append the value of the rho cycle to the list.
        # Increase the number of bits. This meets well with the while loop condition.
        n_bits *= 2

    # If we get no divisor then we return the factorization as it is with a warning that a
    # different rho function may be able to complete the factorization.
    print('Warning unable to decompose value {} using rho pollard factorization.'
          'Try a different rho function'.format(factorization.reduced_value))
    return factorization
