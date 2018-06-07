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
    Dixon
    Garbell_quadratic
    Garbell_cos_nombres
No distinction made between minuscules and capital letters. Default is Garbell_Eratostenes.""",
                        default='Garbell_Eratostenes.')
    args = parser.parse_args()

    # We call the 'factoritza' function that will actually do all the factorizing work.
    factorized = factoritza(args.number, algorithm=args.algorithm)

    # We print the solution.
    print(factorized)









