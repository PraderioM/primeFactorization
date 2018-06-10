import argparse

import factorization.__init__


def get_primality_tests(primality_tests, k=10):
    primality_tests_functions = []
    for test in primality_tests:
        test = test.lower()
        if test == 'trial_division':
            primality_tests_functions.append(factorization.primality_tests.trial_division)
        elif test == 'pseudoprime_primality_test':
            primality_tests_functions.append(lambda n: factorization.primality_tests.pseudoprime_primality_test(n, k))
        elif test == 'solovay_strassen':
            primality_tests_functions.append(lambda n: factorization.primality_tests.solovay_strassen_primality_test(n,
                                                                                                                     k))
        elif test == 'miller_rabin':
            primality_tests_functions.append(lambda n: factorization.primality_tests.miller_rabin_primality_test(n, k))
        else:
            error_msg = 'Unknown primality test {}.'.format(test)
            raise ValueError(error_msg)

    return primality_tests_functions


def factoritza(n, algorithm='Garbell_cos_nombres',
               primality_tests=None, k=10):
    # if the input is 0, 1 or -1 the factorization is trivial and we return the result.
    if primality_tests is None:
        primality_tests = ['miller_rabin']
    if n in [0, 1, -1]:
        return factorization.utils.Factorization(n)

    primality_tests = get_primality_tests(primality_tests, k=k)

    # look for the selected algorithm within all possible algorithms and execute it.
    algorithm = algorithm.lower()
    if algorithm == 'garbell_eratostenes':
        return factorization.garbell_eratostenes.factorize(n)
    elif algorithm == 'rho_pollard':
        return factorization.rho_pollard.factorize(n, primality_tests=primality_tests)
    elif algorithm == 'p-1_pollard':
        return factorization.p1_pollard.factorize(n, primality_tests=primality_tests)
    elif algorithm == 'dixon':
        return factorization.dixon.factorize(n, primality_tests=primality_tests)
    elif algorithm == 'garbell_quadratic':
        return factorization.garbell_quadratic.factorize(n, primality_tests=primality_tests)
    else:
        # if there was no algorithm math we raise a value error
        error_msg = 'Unrecognized factorization algorithm {}.'.format(algorithm)
        raise ValueError(error_msg)


if __name__ == '__main__':
    # Create a parser in order to be able to call this script from a command line giving it some arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--number', type=int, help='Integer to be factorized.')
    parser.add_argument('--algorithm',
                        help="""Factorization algorithm to be used.
It can be one of the following:
    Garbell_Eratostenes
    Rho_Pollard
    p-1_Pollard
    Dixon
    Garbell_quadratic
No distinction made between minuscules and capital letters. Default is Dixon.""",
                        default='dixon')
    parser.add_argument('--primality-test', nargs='+',
                        help="""Primality tests to be used.
It can be one or more of the following:
    trial_division
    pseudoprime_primality_test
    solovay_strassen
    miller_rabin
No distinction made between minuscules and capital letters. Default is miller_rabin.""",
                        default='miller_rabin')
    parser.add_argument('-k', '--rep', type=int, default=10,
                        help='Number of times a primality test should be applied.')
    args = parser.parse_args()

    # We call the 'factoritza' function that will actually do all the factorizing work.
    factorized = factoritza(args.number, algorithm=args.algorithm,
                            primality_tests=[args.primality_test], k=args.rep)

    # We print the solution.
    print(factorized)
