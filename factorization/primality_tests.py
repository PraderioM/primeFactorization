from math import sqrt
from random import randint

import factorization.utils as utils


def is_decomposed(factorization,
                  primality_tests=None):
    """
    Runs a specified list of primality test on the reduced value of a factorization object and uses them to decide if
    the factorization is complete.
    :param factorization: Factorization object that we want to decide if is factorized or not.
    :param primality_tests: List of primality tests to be run to decide if the factorization is completed.
    :return: True if the reduced_value of the factorization object passes all the primality tests. False otherwise.
    """
    # Initialize primality tests.
    if primality_tests is None:
        primality_tests = [miller_rabin_primality_test]

    # Check if the reduced value of the factorization is a prime number by running all primality tests specified.
    reduced_value = factorization.reduced_value
    for primality_test in primality_tests:
        # If the reduced value does not pass one of the primality tests then it is not prime.
        if not primality_test(reduced_value):
            return False
    # If all primality test were passed we can assume the factorization object is factorized.
    return True


def trial_division(n):
    """
    This function tries to divide an integer n by 2, 3, and all integers congruent to 1 or -1 modulo 6
    (candidates to be integers) lower than sqrt(n). If the none of these number is a divisor
    we conclude that n is a prime number and return true.
    :param n: An integer.
    :return: True if n is prime False otherwise.
    """
    n = abs(n)
    # Try for 2 and 3.
    if (n % 2) == 0 or (n % 3) == 0:
        return False

    # Try for integers lower than sqrt(n).
    i = 1
    max_value = sqrt(n) + 1
    # All primes other than 2 and 3 are congruent to 1 or -1 modulo 6.
    while 6 * i - 1 < max_value:
        if (n % (6 * i - 1)) == 0:
            return False
        if (n % (6 * i + 1)) == 0:
            return False

    return True


def is_pseudoprime_to_base_b(n, b=2):
    """
    Checks if a given integer is pseudoprime to a given base or has common factors with that base.
    :param n: an integer.
    :param b: an integer.
    :return: a boolean indicating if n is pseudoprime in the base b.
    """
    if utils.euclidean_algorithm_m_c_d(n, b) != 1:
        return False
    # get the order of b in Z/(n)
    order = utils.get_order(b, n)
    # return True if  the order divides n - 1 and else return False
    return ((n - 1) % order) == 0


def pseudoprime_random_base(n):
    """
    Checks if a given integer is pseudoprime to a random base.
    :param n: an integer.
    :return: a boolean indicating if n is pseudoprime in a random base.
    """
    n = abs(n)
    if n <= 3:  # lower or equal than 3 is prime or 1.
        return True
    b = randint(2, n-1)
    return is_pseudoprime_to_base_b(n, b)


def pseudoprime_primality_test(n, k=10):
    """
    Check if an integer n is a pseudoprime for k random bases.
    :param n: Integer to check primality on.
    :param k: Number of times to check if n if a pseudoprime in a random base.
    :return: True if n passes the pseudoprime primality test.
    """
    for i in range(k):
        if not pseudoprime_random_base(n):
            return False
    return True


def is_euler_pseudoprime_to_base_b(n, b=2):
    """
    Check if a given number is even or is an euler pseudoprime to a given base
    :param n: an integer.
    :param b: an integer.
    :return: a boolean indicating if n is an euler pseudoprime in the base b.
    """
    # If n is even and not 2 then it is not prime.
    if (n % 2) == 0:
        if n == 2:
            return True
        else:
            return False

    # get Jacobi symbol.
    jb_symbol = utils.get_jacobi_symbol(b, n)
    # get b^((n-1)/2) mod n.
    power_mod_n = utils.get_power_mod_n(b, (n-1)/2, n)
    # return True if the two previous numbers are equal modulo n else False.
    return ((jb_symbol - power_mod_n) % n) == 0


def euler_pseudoprime_random_base(n):
    """
    Checks if a given integer is an euler pseudoprime to a random base.
    :param n: an integer.
    :return: a boolean indicating if n is an euler pseudoprime in a random base.
    """
    n = abs(n)
    if n <= 3:   # lower or equal than 3 is prime or 1.
        return True
    b = randint(2, n-1)
    return is_euler_pseudoprime_to_base_b(n, b)


def solovay_strassen_primality_test(n, k=10):
    """
    Check if an integer n is an Euler pseudoprime for k random bases.
    :param n: Integer to check primality on.
    :param k: Number of times to check if n if an euler pseudoprime in a random base.
    :return: True if n passes the Solovay-Strassen primality test else False.
    """
    for i in range(k):
        if not euler_pseudoprime_random_base(n):
            return False
    return True


def is_strong_pseudoprime_to_base_b(n, b=2):
    """
    Check if a given number is even or is an strong pseudoprime to a given base
    :param n: an integer.
    :param b: an integer.
    :return: a boolean indicating if n is an strong pseudoprime in the base b.
    """
    # If n is even it is not prime unless it is 2.
    if (n % 2) == 0:
        if n == 2:
            return True
        else:
            return False

    # we can write n-1 = 2^s*t with t odd.
    t = n - 1
    s = 0
    while (t % 2) == 0:
        s += 1
        t /= 2
    # get b^t
    b = utils.get_power_mod_n(b, t, n)
    if b == 1:
        return True
    for i in range(s):
        # if b is zero then it has no inverse on Z/(n) and, therefore, n is not prime.
        if b == 0:
            return False
        # if b is -1 mod n then it n is a strong pseudoprime in base b
        if b + 1 == n:
            return True
        # if we reach 1 before reaching -1 then n is not a strong pseudoprime in base b.
        if b == 1:
            return False
        # we compute the next square power of b modulo n (b^(2^{i+1}t)).
        b = (b * b) % n

    # This code should not be reached but if -1 was not reached
    # up until now then n is not a strong pseudoprime in base b.
    return False


def strong_pseudoprime_random_base(n):
    """
    Checks if a given integer is a strong pseudoprime to a random base.
    :param n: an integer.
    :return: a boolean indicating if n is a strong pseudoprime in a random base.
    """
    n = abs(n)
    if n <= 3:  # lower or equal than 3 is prime or 1.
        return True
    b = randint(2, n-1)
    return is_strong_pseudoprime_to_base_b(n, b)


def miller_rabin_primality_test(n, k=10):
    """
    Check if an integer n is a strong pseudoprime for k random bases.
    :param n: Integer to check primality on.
    :param k: Number of times to check if n if a strong pseudoprime in a random base.
    :return: True if n passes the Solovay-Strassen primality test else False.
    """
    for i in range(k):
        if not strong_pseudoprime_random_base(n):
            return False
    return True
