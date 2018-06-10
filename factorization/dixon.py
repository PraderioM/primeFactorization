from math import sqrt
from random import randint

import utils
import primality_tests as primality


def get_bi(n, factor_base, max_iters=100):
    """
    Gets a set of
    :param n: A positive integer.
    :param factor_base: A list of prime numbers.
    :param max_iters: Maximum number of iterations to be performed to find the bi values.
    :return: A list of tuples of type (int, list[int]).
            The integers are such that their squares are smooth respect to the factor base.
            The lists are vector representations of that integers square least absolute residue modulo n.
    """
    # initialize bi list.
    bi_list = []
    # initialize vector representations of those bi
    bi_vectors = []
    for i in range(max_iters):
        # look randomly for possible candidates.
        bi_candidate = randint(int(sqrt(n)) - 1, n - 1)
        # Check if the square of the candidate is smooth in the factor base.
        square_residue = utils.least_absolute_residue(bi_candidate ** 2, n)
        if utils.is_smooth(square_residue, factor_base):
            bi_list.append(bi_candidate)
            bi_vectors.append(utils.get_vector_representation(square_residue, factor_base))
            stop, kernel = utils.is_linear_dependent(bi_vectors)
            if stop:
                return [(bi, vector_representation) for bi, vector_representation, index
                        in zip(bi_list, bi_vectors, kernel) if index == 1]


def recursive_step(divisor, factorization, primality_tests, factor_base, max_iters):
    # Get the factorization of the mcd recursively.
    mcd_factorization = factorize(divisor,
                                  primality_tests=primality_tests,
                                  factor_base=factor_base,
                                  max_iters=max_iters)
    # append it to the factorization of n.
    factorization = utils.add_divisor_factorization(factorization, mcd_factorization)

    # we do the same with the current reduced value of the factorization.
    reduced_value_factorization = factorize(factorization.reduced_value,
                                            primality_tests=primality_tests,
                                            factor_base=factor_base,
                                            max_iters=max_iters)
    factorization = utils.add_divisor_factorization(factorization, reduced_value_factorization)

    # We now break the loop and return the factorized value.
    return factorization


def factorize(n, primality_tests=None, factor_base=None, bound=100,
              max_iters=100):
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

    # Get a factor base if needed.
    if factor_base is None:
        factor_base = utils.get_lower_primes(bound)

    for i in range(max_iters):
        # get a list of bi such that their vector representations in the factor_base sum 0 in Z/(2).
        bi_vector_list = get_bi(n, factor_base, max_iters=max_iters)
        if bi_vector_list is None:
            continue
        b_list, vector_list = zip(*bi_vector_list)
        b_list = list(b_list)
        vector_list = list(vector_list)

        # multiply all bi.
        b = utils.get_prod_mod_n(b_list, n)
        b = utils.least_absolute_residue(b, n)  # get the result's least absolute residue.

        # Add together all indexes of the vectors in vector list and divide the result by 2.
        vector_length = len(vector_list[0])
        c_in_factor_base = [sum([vector[i] for vector in vector_list]) / 2 for i in range(vector_length)]
        # Get the element corresponding to the obtained vector in the factor base.
        c = utils.get_element_from_vector(c_in_factor_base, factor_base, n)
        c = utils.least_absolute_residue(c, n)  # get its least absolute residue.
        # This numbers square should be congruent to b modulo n.

        # Check that if b + c has a common non trivial factor with n.
        mcd = utils.euclidean_algorithm_m_c_d(b+c, n)
        if mcd not in [1, n]:
            # We have found a divisor of n by sheer luck.
            return recursive_step(mcd, factorization, primality_tests, factor_base, max_iters)
        # If not check if b - c has a common non trivial factor with n.
        mcd = utils.euclidean_algorithm_m_c_d(b-c, n)
        if mcd not in [1, n]:
            # We have found a divisor of n by sheer luck.
            return recursive_step(mcd, factorization, primality_tests, factor_base, max_iters)
        # If this is not the case either we try again.

    # If we get no divisor then we return the factorization as it is with a warning that a
    # different rho function may be able to complete the factorization.
    print('Warning unable to decompose value {} using dixon factorization.'.format(factorization.reduced_value))
    return factorization
