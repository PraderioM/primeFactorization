import argparse
from datetime import datetime

from factorization.utils import Factorization, get_lower_primes, add_divisor_factorization
from factorization import primality_tests, garbell_eratostenes, rho_pollard, p1_pollard, dixon, garbell_quadratic


def get_primality_tests(tests, k=10):
    """

    :param tests: Names of primality tests to get.
    :param k: number of times each primality test should be tried.
    :return: a list of primality tests.
    """
    primality_tests_functions = []
    for test in tests:
        test = test.lower()
        if test == 'trial_division':
            primality_tests_functions.append(primality_tests.trial_division)
        elif test == 'pseudoprime_primality_test':
            primality_tests_functions.append(lambda n: primality_tests.pseudoprime_primality_test(n, k))
        elif test == 'solovay_strassen':
            primality_tests_functions.append(lambda n: primality_tests.solovay_strassen_primality_test(n, k))
        elif test == 'miller_rabin':
            primality_tests_functions.append(lambda n: primality_tests.miller_rabin_primality_test(n, k))
        else:
            error_msg = 'Unknown primality test {}.'.format(test)
            raise ValueError(error_msg)

    return primality_tests_functions


def find_divisor(n, algorithms):
    n = abs(n)
    # Try to find a divisor using all selected algorithms
    for algorithm in algorithms:
        divisor = algorithm(n)

        if divisor not in [1, n]:
            return divisor

    # If we where unable to find a non trivial divisor we return a trivial one.
    return n


def recursive_step(divisor, factorization, algorithms, tests, max_iters):
    # Get the factorization of the mcd recursively.
    mcd_factorization = recursive_factorization(divisor, algorithms, tests, max_iters)
    # append it to the factorization of n.
    factorization = add_divisor_factorization(factorization, mcd_factorization)

    # we do the same with the current reduced value of the factorization.
    reduced_value_factorization = recursive_factorization(factorization.reduced_value,
                                                          algorithms, tests, max_iters)
    factorization = add_divisor_factorization(factorization, reduced_value_factorization)

    # We now break the loop and return the factorized value.
    return factorization


def recursive_factorization(n, algorithms, tests, max_iters):
    # Initialize factorization object.
    factorization = Factorization(n)
    n = abs(n)  # Set n as a positive number to avoid problems.

    # Check if the Factorization object is already decomposed (if its reduced value is prime).
    # If the reduced value is indeed prime we are done,
    # add the reduced value to the factorization and return the factorization.
    if primality_tests.is_decomposed(factorization, tests):
        factorization.add_factor(factorization.reduced_value)
        return factorization

    for i in range(max_iters):
        divisor = find_divisor(n, algorithms)
        # If we are able to find a non trivial divisor we apply the algorithm recursively.
        if divisor not in [1, n]:
            return recursive_step(divisor, factorization, algorithms, tests, max_iters)

    # If code reaches this far then it means that it was unable to completely factorize input value.
    # Show a warning message.
    print('Warning unable to decompose value {} using selected algorithms.'.format(factorization.reduced_value))
    return factorization


def get_find_divisor_algorithms(algorithms, prime_list):
    """
    :param algorithms: A list of algorithm names.
    :param prime_list: A list of primes or None.
    :return: A list of functions that take as input a integer value and output a divisor of such integer.
    """
    find_divisor_functions = []
    for algorithm in algorithms:
        if algorithm == 'rho_pollard':
            find_divisor_functions.append(lambda n: rho_pollard.find_divisor(n, max_iters=1))
        elif algorithm == 'p-1_pollard':
            find_divisor_functions.append(lambda n: p1_pollard.find_divisor(n, max_iters=1))
        elif algorithm == 'dixon':
            find_divisor_functions.append(lambda n: dixon.find_divisor(n, factor_base=prime_list,
                                                                       max_iters=1))
        elif algorithm == 'garbell_quadratic':
            find_divisor_functions.append(lambda n: garbell_quadratic.find_divisor(n, small_primes=prime_list,
                                                                                   max_iters=1))
        elif algorithm == 'garbell_eratostenes':
            find_divisor_functions.append(lambda n: garbell_eratostenes.find_divisor(n))
        else:
            # if there was no algorithm math we raise a value error
            error_msg = 'Unrecognized factorization algorithm {}.'.format(algorithm)
            raise ValueError(error_msg)
    return find_divisor_functions


def factoritza(n, algorithms=None,
               tests=None, k=10, max_iters=100, low_primes_bound=5000):
    """
    :param n: number to be factorized.
    :param algorithms: Algorithms to be used for factorization. If None they will be set to
                       ['Garbell_Eratostenes', 'p-1_Pollard', 'Rho_Pollard']
    :param tests: Primality tests to be used to test primality of factors. If None it will be set to
                            ['miller_rabin']
    :param k: Number of times a primality test should be applied.
    :param max_iters: number of times a factorization algorithm should be applied at every factorization step.
    :param low_primes_bound: Bound for lower primes to be removed with Eratostenes Sieve algorithm.
    :return: A Factorize object containing the factorized integer.
    """

    if algorithms is None:
        algorithms = ['Garbell_Eratostenes', 'p-1_Pollard', 'Rho_Pollard']
    if tests is None:
        tests = ['miller_rabin']

    # if the input is 0, 1 or -1 the factorization is trivial and we return the result.
    if n in [0, 1, -1]:
        print('Factorization_start')
        start_time = datetime.now()
        result = Factorization(n)
        end_time = datetime.now()
        print('Factorization completed in {} seconds.'.format((end_time - start_time).total_seconds()))
        return result

    # If dixon or garbell_quadratic algorithm are in algorithm list then we need a prime list.
    # We should look for it else it will be None.
    if 'dixon' in algorithms or 'garbell_quadratic' in algorithms:
        prime_list = get_lower_primes(bound=low_primes_bound)
    else:
        prime_list = None

    # we get functions corresponding to the selected primality tests and algorithms.
    print('Initializing primality tests.')
    tests = get_primality_tests(tests, k=k)
    print('Initializing algorithms.')
    find_divisor_algorithms = get_find_divisor_algorithms(algorithms, prime_list)

    # Factorization start.
    print('factorization_start')
    start_time = datetime.now()

    # Execute first the Garbell_eratostenes algorithm to remove small primes and find a list of small primes.
    # If it is the only selected algorithm we return the factorized number.
    if 'garbell_eratostenes' in algorithms:
        algorithms.remove('garbell_eratostenes')
        if len(algorithms) == 0:
            factorization = garbell_eratostenes.factorize(n, max_value=None)
            end_time = datetime.now()
            print('Factorization completed in {} seconds.'.format((end_time - start_time).total_seconds()))
            return factorization
        else:
            factorization = garbell_eratostenes.factorize(n, max_value=low_primes_bound)
    else:
        # Initialize factorization.
        factorization = Factorization(n)

    # Start factorization for obtained reduced value.
    reduced_value_factorization = recursive_factorization(factorization.reduced_value,
                                                          algorithms=find_divisor_algorithms,
                                                          tests=tests,
                                                          max_iters=max_iters)
    # Add new factors.
    factorization = add_divisor_factorization(factorization, reduced_value_factorization)

    end_time = datetime.now()
    print('Factorization completed in {} seconds.'.format((end_time - start_time).total_seconds()))
    return factorization


if __name__ == '__main__':
    # Create a parser in order to be able to call this script from a command line giving it some arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--number', type=int, help='Integer to be factorized.')
    parser.add_argument('-mi', '--max-iters', type=int, default=100,
                        help='Maximum number of times the factorization algorithm should be tried'
                             'at every factorization step. Default is 100.')
    parser.add_argument('--algorithm', nargs='+',
                        help="Factorization algorithms to be used.\n"
                             "It can be one of the following:\n"
                             "\tGarbell_Eratostenes\n"
                             "\tRho_Pollard\n"
                             "\tp-1_Pollard\n"
                             "\tDixon\n"
                             "\tGarbell_quadratic\n"
                             "No distinction made between minuscules and capital letters.\n"
                             "Default is Garbell_Eratostenes Rho_Pollard p-1_Pollard.",
                        default='Garbell_Eratostenes p-1_Pollard Rho_Pollard')
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
    factorized = factoritza(args.number, algorithms=[args.algorithm],
                            tests=[args.primality_test], k=args.rep,
                            max_iters=args.max_iters)

    # We print the solution.
    print(factorized)
