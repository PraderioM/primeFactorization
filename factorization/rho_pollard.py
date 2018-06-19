from random import randint

import factorization.utils as utils
import factorization.primality_tests as primality


def recursive_step(divisor, factorization,
                   primality_tests,
                   rho_function):
    # Get the factorization of the mcd recursively.
    mcd_factorization = factorize(divisor, primality_tests=primality_tests,
                                  rho_function=rho_function)
    # append it to the factorization of n.
    factorization = utils.add_divisor_factorization(factorization, mcd_factorization)

    # we do the same with the current reduced value of the factorization.
    reduced_value_factorization = factorize(factorization.reduced_value,
                                            primality_tests=primality_tests,
                                            rho_function=rho_function)
    factorization = utils.add_divisor_factorization(factorization, reduced_value_factorization)

    # We now break the loop and return the factorized value.
    return factorization


def factorize(n, primality_tests=None,
              rho_function=None, max_iters=10000):
    """
    Applies the rho pollard factorization method
    :param n: Integer to be factorized.
    :param primality_tests: Primality tests to be used.
    :param rho_function: rho function to be applied.
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

    # Initialize rho_function.
    if rho_function is None:
        # Random non zero integer for rho function.
        aux = randint(1, n - 1)

        # Quadratic rho function.
        def rho_function(x):
            return (x ** 2 + aux) % n

    # Try the algorithm at most max_iters times.
    for i in range(max_iters):
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
                mcd = utils.euclidean_algorithm_m_c_d(n, last_rho_function_image-compare_value)
                # if the mcd is not trivial we continue have found a divisor and we continue the algorithm recursively.
                if mcd not in [1, n]:
                    return recursive_step(mcd, factorization, primality_tests, rho_function)

                # Else we append the value of the rho cycle to the list.
            # Increase the number of bits. This meets well with the while loop condition.
            n_bits *= 2

    # If we get no divisor then we return the factorization as it is with a warning that a
    # different rho function may be able to complete the factorization.
    print('Warning unable to decompose value {} using rho pollard factorization.'.format(factorization.reduced_value))
    print('Try a different rho function')
    return factorization
