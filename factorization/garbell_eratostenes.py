from math import sqrt

import factorization.utils as utils


def factorize(n):
    """
    Applies the sieve of eratostenes and fermat to decompose a number in its prime factors.
    :param n:
    :return:
    """
    factorization = utils.Factorization(n)

    # initialize the sieve.
    max_value = int(sqrt(factorization.reduced_value))+1
    sieve = [i for i in range(2, max_value)]
    # apply_algorithm iteratively
    while len(sieve) != 0:
        # get first element of the sieve (it is a prime number).
        p = sieve.pop(0)
        # If it is possible divide by the prime p.
        while factorization.reduced_value % p == 0:
            factorization.add_factor(p)

        # remove all values of the sieve divisible by p and greater than the square root of the new reduced value.
        max_value = int(sqrt(factorization.reduced_value)) + 1
        sieve = [i for i in sieve if i % p != 0 and i < max_value]

    # now that the sieve is empty it can either mean two things.
    # The first that the reduced value is now a prime number the second that it is 1.
    if factorization.reduced_value != 1:
        factorization.add_factor(factorization.reduced_value)

    return factorization
