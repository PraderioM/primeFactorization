from math import sqrt

import factorization.utils as utils


def factorize(n, max_value=5000, return_found_primes=False):
    """
    Applies the sieve of eratostenes and fermat to decompose a number in its prime factors.
    :param n: An integer to be factorized.
    :param max_value: An integer indicating the maximum value to be used in the sieve for factorizing.
                      If None it will be set to sqrt(n)
    :param return_found_primes: Boolean indicating if the list of primes found in the factorization process
                                should be returned.
    :return: Factorization object with the factorized input integer.
    """
    factorization = utils.Factorization(n)

    # Initialize max value if not already initialized.
    if max_value is None:
        last_is_prime = True
        max_value = int(sqrt(factorization.reduced_value))+1
    else:
        last_is_prime = False
        max_value = min(max_value, int(sqrt(factorization.reduced_value)) + 1)

    # initialize the sieve.
    sieve = [i for i in range(2, max_value)]
    found_primes = []
    # apply_algorithm iteratively.
    while len(sieve) != 0:
        # get first element of the sieve (it is a prime number).
        p = sieve.pop(0)
        found_primes.append(p)
        # If it is possible divide by the prime p.
        while factorization.reduced_value % p == 0:
            factorization.add_factor(p)

        # remove all values of the sieve divisible by p and greater than the square root of the new reduced value.
        max_value = min(int(sqrt(factorization.reduced_value)) + 1, max_value)
        sieve = [i for i in sieve if i % p != 0 and i < max_value]

    # now that the sieve is empty it can either mean two things if max_value was originally None.
    # The first that the reduced value is now a prime number the second that it is 1.
    if last_is_prime:
        factorization.add_factor(factorization.reduced_value)

    # Return all prime numbers found if asked so.
    if return_found_primes:
        return factorization, found_primes
    return factorization


def find_divisor(n):
    """
    Applies the sieve of eratostenes and fermat find a prime divisor for n.
    :param n: An integer to be factorized.
    :return: A prime divisor of n.
    """
    n = abs(n)  # Set n as a positive number to avoid problems.

    # Initialize max value if not already initialized.
    max_value = int(sqrt(n))+1

    # initialize the sieve.
    sieve = [i for i in range(2, max_value)]
    # apply_algorithm iteratively.
    while len(sieve) != 0:
        # get first element of the sieve (it is a prime number).
        p = sieve.pop(0)
        # If p is a divisor of n we return p.
        if n % p == 0:
            return p

        sieve = [i for i in sieve if i % p != 0 and i < max_value]

    # This code should never be reached but just to make sure we return a trivial divisor of n.
    return n
