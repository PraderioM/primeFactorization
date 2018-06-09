import argparse

from factorization.utils import factoritza


def main():
    # Create a parser in order to be able to call this script from a command line giving it some arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--number', type=int, help='Integer to be factorized.')
    parser.add_argument('--algorithm',
                        help="""Factorization algorithm to be used.
It can be one of the following:
    Garbell_Eratostenes
    Rho_Pollard
    p-1_Pollard
    Garbell_quadratic
    Garbell_cos_nombres
No distinction made between minuscules and capital letters. Default is Garbell_cos_nombres.""",
                        default='Garbell_cos_nombres.')
    parser.add_argument('--primality-test', nargs='+',
                        help="""Primality tests to be used.
It can be one or more of the following:
    trial_division
    pseudoprime_primality_test
    solovay_strassen
    miller_rabin
No distinction made between minuscules and capital letters. Default is miller_rabin.""",
                        default='miller_rabin.')
    parser.add_argument('-k', '--rep', type=int, default=10,
                        help='Number of times a primality test should be applied.')
    args = parser.parse_args()

    # We call the 'factoritza' function that will actually do all the factorizing work.
    factorized = factoritza(args.number, algorithm=args.algorithm,
                            primality_tests=args.primality_test, k=args.rep)

    # We print the solution.
    print(factorized)









