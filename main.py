import argparse
from datetime import datetime

import factorization.primality_tests


def get_primality_tests(primality_tests, k=10):
    """

    :param primality_tests:
    :param k:
    :return:
    """
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
               primality_tests=None, k=10, max_iters=10000):
    """
    :param n: number to be factorized.
    :param algorithm: algorithm to be used for factorization.
    :param primality_tests: primality tests to be used to test primality of factors.
    :param k: Number of times a primality test should be applied.
    :param max_iters: number of times a factorization algorithm should be applied at every factorization step.
    :return: A Factorize object containing the factorized integer.
    """
    # if the input is 0, 1 or -1 the factorization is trivial and we return the result.
    if primality_tests is None:
        primality_tests = ['miller_rabin']
    if n in [0, 1, -1]:
        import factorization.utils
        print('Factorization_start')
        start_time = datetime.now()
        result = factorization.utils.Factorization(n)
        end_time = datetime.now()
        print('Factorization completed in {} seconds.'.format((end_time-start_time).total_seconds()))
        return result

    primality_tests = get_primality_tests(primality_tests, k=k)

    # look for the selected algorithm within all possible algorithms and execute it.
    algorithm = algorithm.lower()
    if algorithm == 'garbell_eratostenes':
        import factorization.garbell_eratostenes
        print('factorization_start')
        start_time = datetime.now()
        result = factorization.garbell_eratostenes.factorize(n)
        end_time = datetime.now()
        print('Factorization completed in {} seconds.'.format((end_time-start_time).total_seconds()))
        return result
    elif algorithm == 'rho_pollard':
        import factorization.rho_pollard
        print('factorization_start')
        start_time = datetime.now()
        result = factorization.rho_pollard.factorize(n, primality_tests=primality_tests, max_iters=max_iters)
        end_time = datetime.now()
        print('Factorization completed in {} seconds.'.format((end_time-start_time).total_seconds()))
        return result
    elif algorithm == 'p-1_pollard':
        import factorization.p1_pollard
        print('factorization_start')
        start_time = datetime.now()
        result = factorization.p1_pollard.factorize(n, primality_tests=primality_tests, max_iters=max_iters)
        end_time = datetime.now()
        print('Factorization completed in {} seconds.'.format((end_time-start_time).total_seconds()))
        return result
    elif algorithm == 'dixon':
        import factorization.dixon
        print('factorization_start')
        start_time = datetime.now()
        result = factorization.dixon.factorize(n, primality_tests=primality_tests, max_iters=max_iters)
        end_time = datetime.now()
        print('Factorization completed in {} seconds.'.format((end_time-start_time).total_seconds()))
        return result
    elif algorithm == 'garbell_quadratic':
        import factorization.garbell_quadratic
        print('factorization_start')
        start_time = datetime.now()
        result = factorization.garbell_quadratic.factorize(n, primality_tests=primality_tests, max_iters=max_iters)
        end_time = datetime.now()
        print('Factorization completed in {} seconds.'.format((end_time-start_time).total_seconds()))
        return result
    else:
        # if there was no algorithm math we raise a value error
        error_msg = 'Unrecognized factorization algorithm {}.'.format(algorithm)
        raise ValueError(error_msg)


if __name__ == '__main__':
    # Create a parser in order to be able to call this script from a command line giving it some arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--number', type=int, help='Integer to be factorized.')
    parser.add_argument('-mi', '--max-iters', type=int, default=10000,
                        help='Maximum number of times the factorization algorithm should be tried'
                             'at every factorization step. Default is 10,000.')
    parser.add_argument('--algorithm',
                        help="Factorization algorithm to be used.\n"
                             "It can be one of the following:\n"
                             "\tGarbell_Eratostenes\n"
                             "\tRho_Pollard\n"
                             "\tp-1_Pollard\n"
                             "\tDixon\n"
                             "\tGarbell_quadratic\n"
                             "No distinction made between minuscules and capital letters. Default is p-1_Pollard.",
                        default='p-1_Pollard')
    parser.add_argument('--primality-test', nargs='+',
                        help="Primality tests to be used.\n"
                             "It can be one or more of the following:\n"
                             "\ttrial_division\n"
                             "\tpseudoprime_primality_test\n"
                             "\tsolovay_strassen\n"
                             "\tmiller_rabin\n"
                             "No distinction made between minuscules and capital letters. Default is miller_rabin.",
                        default='miller_rabin')
    parser.add_argument('-k', '--rep', type=int, default=10,
                        help='Number of times a primality test should be applied.')
    args = parser.parse_args()

    # We call the 'factoritza' function that will actually do all the factorizing work.
    factorized = factoritza(args.number, algorithm=args.algorithm,
                            primality_tests=[args.primality_test], k=args.rep, max_iters=args.max_iters)

    # We print the solution.
    print(factorized)
