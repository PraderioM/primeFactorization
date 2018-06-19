from math import sqrt, exp, log
from random import randint

import utils
import primality_tests as primality


def solve_x2_equal_a_mod_2k(a, k):
    if k in [1, 2, 3]:
        return 1
    else:
        # initialize for k=3.
        xk = 1
        k2 = 8
        for i in range(k - 3):
            aux = 0 if ((xk ** 2 - a) / k2) % 2 == 0 else 1
            # Get solution for k+1.
            xk = xk + aux * k2 / 2
            k2 *= 2
        return xk


def tonelli_shanks(a, p):
    """
    It is supposed that a solution exists and we do not have a closed formula for it.
    Square root modulo prime number
    Solve the equation
        x^2 = a mod p
    and return list of x solution
    http://en.wikipedia.org/wiki/Tonelli-Shanks_algorithm
    """
    a %= p

    # Simple case
    if a == 0:
        return 0
    if p == 2:
        return a

    # Factor p-1 on the form q * 2^s (with q odd).
    q, s = p - 1, 0
    while q % 2 == 0:
        s += 1
        q /= 2

    # Select a z which is a quadratic non residue modulo p.
    z = 1
    while utils.get_jacobi_symbol(z, p) != -1:
        z += 1
    c = utils.get_power_mod_n(z, q, p)

    # Search for a solution
    x = utils.get_power_mod_n(a, (q + 1) / 2, p)
    t = utils.get_power_mod_n(a, q, p)
    m = s
    while t != 1:
        # Find the lowest i such that t^(2^i) = 1
        i, e = 0, 2
        for i in range(1, m):
            if utils.get_power_mod_n(t, e, p) == 1:
                break
            e *= 2

        # Update next value to iterate
        b = utils.get_power_mod_n(c, 2 ** (m - i - 1), p)
        x = (x * b) % p
        t = (t * b * b) % p
        c = (b * b) % p
        m = i

    return x


def solve_x2_equal_a_mod_p_odd(a, p):
    # check if a solution exists.
    if utils.get_power_mod_n(a, (p - 1) / 2, p) != 1:
        return None
    # If exists we get it in two different ways.
    # If p % 3 mod 4 we have the following closed formula.
    if (p - 3) % 4 == 0:
        return utils.get_power_mod_n(a, (p + 1) / 4, p)
    # If p % 1 mod 4 we apply Tonelli-Shanks algorithm.
    return tonelli_shanks(a, p)


def recursive_solve_x2_equal_a_mod_pk1_odd(a, pk, solk):
    yk = utils.get_power_mod_n(2 * solk, pk - 1, pk)
    if yk == 0:  # if yk has no inverse there is no solution.
        return None
    else:
        # get solution to x**2=a mod p^(k+1)
        return solk - (solk ** 2 - a) * yk


def solve_x2_equal_a_mod_pk(a, p, k):
    if p % 2 == 0:
        return solve_x2_equal_a_mod_2k(a, k)
    else:
        pk = p
        # get solution to x**2=a mod p^(k+1)
        xk = solve_x2_equal_a_mod_p_odd(a, p)
        for i in range(k - 1):
            # get (2*xk)^(-1)
            xk = recursive_solve_x2_equal_a_mod_pk1_odd(a, pk, xk)
            if xk is None:
                return None
            pk *= p

        return xk


def get_bi(n, factor_base, P):
    """
    Gets a set of
    :param n: A positive integer.
    :param factor_base: A list of prime numbers.
    :param P: P threshold of the quadratic sieve algorithm.
    :return: A list of tuples of type (int, list[int]).
            The integers are such that their squares are smooth respect to the factor base.
            The lists are vector representations of that integers square least absolute residue modulo n.
    """
    # Get the A value of the quadratic sieve algorithm.
    A = randint(P, P ** 2)

    # initialize bi list.
    int_sqrt_n = int(sqrt(n))
    bi_candidates = [t for t in range(int_sqrt_n + 1, int_sqrt_n + A + 1)]

    # For every number in the factor base and every element t in the bi_list we get
    # the maximum beta such that t**2-n=n mod p^beta
    vectors = {bi: {'square': bi**2-n,
                    'vector': []} for bi in bi_candidates}
    for p in factor_base:
        for bi in bi_candidates:
            vectors[bi]['vector'].append(0)

        multiple_bi = [bi for bi in bi_candidates if vectors[bi]['square'] % p == 0]
        while len(multiple_bi) != 0:
            # Add one to the vector representation.
            for bi in multiple_bi:
                vectors[bi]['square'] /= p
                vectors[bi]['vector'][-1] += 1
            # Reduce the list of multiples of p.
            multiple_bi = [bi for bi in multiple_bi if vectors[bi]['square'] % p == 0]

    # Get vector the bi candidates whose square is actually smooth with respect to the factor base.
    bi_list = []
    bi_vectors = []
    for bi in bi_candidates:
        if vectors[bi]['square'] == 1:
            bi_list.append(bi)
            bi_vectors.append(vectors[bi]['vector'])

    is_dependent, kernel = utils.is_linear_dependent(bi_vectors)
    if is_dependent:
        return [(bi, vector_representation) for bi, vector_representation, index
                in zip(bi_list, bi_vectors, kernel) if index == 1]
    else:
        return None


def recursive_step(divisor, factorization, primality_tests, small_primes, max_iters):
    # Get the factorization of the mcd recursively.
    mcd_factorization = factorize(divisor,
                                  primality_tests=primality_tests,
                                  small_primes=small_primes,
                                  max_iters=max_iters)
    # append it to the factorization of n.
    factorization = utils.add_divisor_factorization(factorization, mcd_factorization)

    # we do the same with the current reduced value of the factorization.
    reduced_value_factorization = factorize(factorization.reduced_value,
                                            primality_tests=primality_tests,
                                            small_primes=small_primes,
                                            max_iters=max_iters)
    factorization = utils.add_divisor_factorization(factorization, reduced_value_factorization)

    # We now break the loop and return the factorized value.
    return factorization


def factorize(n, primality_tests=None, small_primes=None,
              max_iters=10000):
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

    # Get the P value of the quadratic sieve algorithm.
    P = int(exp(sqrt(log(n) * log(log(n)))))

    # Get small primes lower than P.
    if small_primes is None:
        small_primes = utils.get_lower_primes(P)
    else:
        small_primes = [p for p in small_primes if p <= P]

    # Remove from small_primes all odd primes such that the legendre number (n/p) is not 1.
    factor_base = [2] + [p for p in small_primes[1:] if utils.get_jacobi_symbol(n, p) == 1]

    # Make sure none of these primes divides n.
    for factor in factor_base:
        if n % factor == 0:
            return recursive_step(factor, factorization, primality_tests, small_primes, max_iters)

    for i in range(max_iters):
        # get a list of bi such that their vector representations in the factor_base sum 0 in Z/(2).
        bi_vector_list = get_bi(n, factor_base, P)
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
        mcd = utils.euclidean_algorithm_m_c_d(b + c, n)
        if mcd not in [1, n]:
            # We have found a divisor of n by sheer luck.
            return recursive_step(mcd, factorization, primality_tests, small_primes, max_iters)
        # If not check if b - c has a common non trivial factor with n.
        mcd = utils.euclidean_algorithm_m_c_d(b - c, n)
        if mcd not in [1, n]:
            # We have found a divisor of n by sheer luck.
            return recursive_step(mcd, factorization, primality_tests, small_primes, max_iters)
        # If this is not the case either we try again.

    # If we get no divisor then we return the factorization as it is with a warning that a
    # different rho function may be able to complete the factorization.
    print('Warning unable to decompose value {} using quadratic sieve.'.format(factorization.reduced_value))
    return factorization
