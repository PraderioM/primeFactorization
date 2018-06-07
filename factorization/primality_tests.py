from math import sqrt
from random import randint

from factorization.utils import euclidean_algorithm_m_c_d


def get_order(a: int, n: int) -> int:
    """
    Computes the order of the integer a modulo in Z/(n).
    :param a: an integer.
    :param n: an integer.
    :return: order of a in Z/(n) or -1 if a has no inverse in Z/(n)
    """
    a = abs(a)
    n = abs(n)
    order = 1
    prod = a % n
    while prod != 1:
        # if prod is 0 then a has no inverse in Z/(n) and we stop the algorithm.
        if prod == 0:
            return -1
        # compute always the power modulo n in order to avoid overflow.
        prod = prod * a % n
        order += 1

    return 1


def get_jacobi_symbol(a: int, n: int) -> int:
    """
    :param a: a positive integer (we don't know how to compute fast the Jacobi symbol for negative integers).
    :param n: a positive odd integer.
    :return: the Jacobi symbol (a/n)
    """
    if a < 0:
        raise ValueError('a should be non negative.')
    if n % 2 == 0:
        raise ValueError('n should be odd but it is not.')

    # get the residue of a modulo n.
    a = a % n
    # if it is zero then the Jacobi symbol is 0.
    if a == 0:
        return 0
    # if it is 1 then the Jacobi symbol is 1.
    if a == 1:
        return 1

    # Extract all powers of 2 multiplying a and compute the Jacobi symbol for those.
    prod = 1
    ratio = 1 if n - 1 % 8 in [0, 6] else -1  # (2/n)=(-1)^((n^2-1)/8)
    while a % 2 == 0:
        a /= 2
        prod *= ratio

    # Now that a is odd we can apply the quadratic reciprocity.
    prod = prod if (a - 1) % 4 == 0 or (n - 1) % 4 == 0 else -prod  # we combine it with the previous results.
    return prod * get_jacobi_symbol(n, a)  # apply the algorithm recursively.


def get_power_mod_n(base, power, mod):
    # We apply an optimized algorithm to compute integer powers in logarithmic time.
    prod = 1
    ratio = base % mod  # do this to avoid overflow
    # while we write the power in binary as power=i_ni_{n-1}...i_1 where i_j are either 0 or 1
    # we compute the product as b^(i_0)(b^2)^(i_1)...(b^{n+1})^(i_n).
    while power != 0:
        if power % 2 == 1:
            prod *= ratio
            prod = prod % mod  # do this to avoid overflow.
            power -= 1
        # increase the ration and divide the power by 2 to get the next binary digit.
        ratio *= base
        ratio = ratio % mod  # do this to avoid overflow.
        power /= 2
    return prod


def trial_division(n: int) -> bool:
    """
    This function tries to divide an integer n by 2, 3, and all integers congruent to 1 or -1 modulo 6
    (candidates to be integers) lower than sqrt(n). If the none of these number is a divisor
    we conclude that n is a prime number and return true.
    :param n: An integer.
    :return: True if n is prime False otherwise.
    """
    n = abs(n)
    # Try for 2 and 3.
    if n % 2 == 0 or n % 3 == 0:
        return False

    # Try for integers lower than sqrt(n).
    i = 1
    max_value = sqrt(n) + 1
    while 6 * i - 1 < max_value:
        if n % (6 * i - 1) == 0:
            return False
        if n % (6 * i + 1) == 0:
            return False

    return True


def is_pseudoprime_to_base_b(n: int, b: int=2) -> bool:
    """
    Checks if a given integer is pseudoprime to a given base or has common factors with that base.
    :param n: an integer.
    :param b: an integer.
    :return: a boolean indicating if n is pseudoprime in the base b.
    """
    if euclidean_algorithm_m_c_d(n, b) == 1:
        return False
    # get the order of b in Z/(n)
    order = get_order(b, n)
    # return True if  the order divides n - 1 and else return False
    return (n - 1) % order == 0


def pseudoprime_random_base(n: int) -> bool:
    """
    Checks if a given integer is pseudoprime to a random base.
    :param n: an integer.
    :return: a boolean indicating if n is pseudoprime in a random base.
    """
    n = abs(n)
    if n <= 1:
        return True
    b = randint(0, n-1)
    return is_pseudoprime_to_base_b(n, b)


def is_euler_pseudoprime_to_base_b(n: int, b: int=2) -> bool:
    """
    Check if a given number is even or is an euler pseudoprime to a given base
    :param n: an integer.
    :param b: an integer.
    :return: a boolean indicating if n is an euler pseudoprime in the base b.
    """
    if n % 2 == 0:
        return False

    # get Jacobi symbol.
    jb_symbol = get_jacobi_symbol(b, n)
    # get b^((n-1)/2) mod n.
    power_mod_n = get_power_mod_n(b, (n-1)/2, n)
    # return True if the two previous numbers are equal modulo n else False.
    return (jb_symbol - power_mod_n) % n == 0


def euler_pseudoprime_random_base(n: int) -> bool:
    """
    Checks if a given integer is an euler pseudoprime to a random base.
    :param n: an integer.
    :return: a boolean indicating if n is an euler pseudoprime in a random base.
    """
    n = abs(n)
    if n <= 1:
        return True
    b = randint(0, n-1)
    return is_euler_pseudoprime_to_base_b(n, b)
