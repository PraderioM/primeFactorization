from typing import Dict, Optional, Any, List, Callable

import factorization.garbell_eratostenes
import factorization.rho_pollard
import factorization.p1_pollard
import factorization.garbell_quadratic
import factorization.garbell_cos_nombres
import factorization.primality_tests


def _is_type_error(value, o_type, type_name):
    if not isinstance(value, o_type):
        error_msg = 'Expected {} but got {}.'.format(type_name, type(value))
        raise TypeError(error_msg)


class Factorization(object):
    """
    This class's purpose will be that of storing the factorization of a given integer
    as well as modifying making it easier to print the factorized number.
    """
    def __init__(self, value, factors=None):
        """
        Initializer.
        :param value: Integer to be factorized.
        :param factors: Dictionary having as keys the factors of the prime decomposition and as values its powers.
        """
        if factors is None:
            factors = {}
        
        # check that the inputs are correct and store them.
        # value.
        self._is_valid_value(value)
        self._value = abs(value)
        # factors.
        self._is_valid_factors(factors)
        self._factors = factors

        # Store sign of input value, check if value is zero and store reduced value.
        self._unit = 1 if value >= 0 else -1
        self._is_zero = value == 0
        self._reduced_value = self._value / self.multiply_factors()

    def __str__(self):
        """
        :return: A string showing the decomposition of value as a product of primes.
        """
        # if value is zero, 1 or -1 we return the trivial factorization.
        if self._is_zero:
            return 'The number 0 factorizes as 0.'
        elif self._value == 1:
            return 'The number {0} factorizes as {0}.'.format(self._unit)

        # get the original input number.
        number = self._value * self._unit

        aux_str = 'in prime factors ' if self._reduced_value == 1 else ''
        out_str = 'The integer {} can be decomposed {}as:\n'.format(number, aux_str)

        if self._unit == -1:
            out_str += '-1*'

        # store factorization in a list and sort it for better visualization
        factors = [(factor, self._factors[factor]) for factor in self._factors]
        factors.sort(key=lambda factor: factor[0])

        # actual construction of the output string.
        for factor, power in factors:
            # Add the factor.
            out_str += str(factor)
            # Add the power if it is non trivial.
            if power != 1:
                out_str += '^{}'.format(power)
            # Add the '*' symbol for concatenating with the next factor.
            out_str += '*'

        # remove last '*' symbol.
        out_str = out_str[:-1]
        return out_str

    # region properties
    # getter for unit.
    @property
    def unit(self):
        return self._unit

    # getter for is_zero.
    @property
    def is_zero(self):
        return self._is_zero

    # getter for reduced_value.
    @property
    def reduced_value(self):
        return self._reduced_value

    # getter and setter for value.
    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value):
        self._is_valid_value(value)
        self._value = abs(value)

        # Store sign of input value and check if it is zero.
        self._unit = 1 if value >= 0 else -1
        self._is_zero = value == 0

    # getter and setter for factors.
    @property
    def factors(self):
        return self._factors

    @factors.setter
    def factors(self, factors):
        self._is_valid_factors(factors)
        self._factors = factors
        self._reduced_value = self._value / self.multiply_factors()

    # endregion

    # region input_validations
    @staticmethod
    def _is_valid_value(value):
        """
        Check if value is a valid input factor.
        :param value: integer.
        :return: True if value is an integer else raises an error.
        """
        _is_type_error(value, int, 'integer')
        # If an error was not raised the input is correct.
        return True

    def _is_valid_factors(self, factors):
        """
        Check if factors is a valid input factor.
        :param factors: Dictionary containing integers as key and their powers as values.
        :return: True if factors are valid input else raises an error.
        """
        # Check if input is a dictionary.
        _is_type_error(factors, dict, 'dictionary')
        # check if keys are integers and values are positive integers.
        for key in factors:
            # both key and value should be positive integers
            value = factors[key]
            key -= 1
            for x in [key, value]:
                _is_type_error(x, int, 'integer')
                if x <= 0:
                    error_msg = 'Key should be greater than 1 and value should be positive but got {}.'.format(x)
                    raise ValueError(error_msg)

        # check if factors divide value.
        if self._value % self.multiply_factors(factors) != 0:
            raise ValueError('Factors do not divide value.')

        # If an error was not raised the input is correct.
        return True

    # endregion
    
    def add_factor(self, value):
        """
        Add a factor to the decomposition and update the reduced value.
        :param value: factor to be added to the decomposition.
        """
        # check if factor is an integer. If not we raise an error.
        _is_type_error(value, int, 'integer')
        # check if factor is positive. If not we raise an error.
        if value <= 0:
            error_msg = 'New factor should be greater than 0 but got {}.'.format(value)
            raise ValueError(error_msg)
        # check if factor is already in the decomposition add if it is not and add one to its power if it is.
        if self._reduced_value % value != 0:
            raise ValueError('Factor does not divide value.')

        # If the factor is 1 we do nothing.
        if value == 1:
            return

        # add the factor.
        if value in self.factors:
            self.factors[value] += 1
        else:
            self.factors[value] = 1

        # update the reduced value.
        self._reduced_value /= value

    def add_factors(self, factors):
        """
        Adds a list of factors to the factorization.
        :param factors: List of factors to be added.
        """
        for factor in factors:
            self.add_factor(factor)

    def multiply_factors(self, factors=None):
        """
        Multiplies object factors and returns the result.
        :param factors: Dictionary containing integers as key and their powers as values.
        :return: the product of the integers of the dictionary raised to the corresponding power.
        """
        # if the factors parameter is none then the object factors will be used.
        if not factors:
            factors = self._factors

        # multiply all factors to the given power and return the result.
        prod = 1
        for factor in factors:
            prod *= pow(factor, factors[factor])

        return prod


def add_divisor_factorization(factorization1, factorization2):
    """
    Improves the factorization of a number via de factorization of a divisor.
    :param factorization1: Factorization objects whose factorization should be improved.
    :param factorization2: Factorization of a divisor of the value of the first factorization object.
    :return: Improved factorization.
    """
    # First check if the second factorization value actually divides the first factorization reduced value.
    if factorization2.value % factorization1.reduced_value != 0:
        raise ValueError('First second factorization value does not divide first factorization reduced value')

    # Extract a list of integers from the dictionary of factors.
    power_factors = factorization2.factors
    factors = []
    for factor in power_factors:
        factors.extend([factor]*power_factors[factor])

    # Add every factor of the second factorization factors.
    factorization1.add_factors(factors)

    # We also add the reduced value of the factorization just in case the factorization was not completed.
    factorization1.add_factor(factorization2.reduced_value)

    return factorization1


def euclidean_algorithm_m_c_d(n, m):
    """
    Applies the euclidean algorithm to find the m.c.d. between two integers.
    :param n: an integer.
    :param m: an integer.
    :return: m.c.d. between n and m.
    """
    # make sure numbers are positive.
    n = abs(n)
    m = abs(m)
    if m == 0:
        return n
    # apply the algorithm until we find a residual equal to 0.
    while n % m != 0:
        r = n % m  # r is the residual of the division between n and m.
        n = m
        m = r

    # once the algorithm is completed m is the m.c.d
    return m


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
        return Factorization(n)

    primality_tests = get_primality_tests(primality_tests, k=k)

    # look for the selected algorithm within all possible algorithms and execute it.
    algorithm = algorithm.lower()
    if algorithm == 'garbell_eratostenes':
        return factorization.garbell_eratostenes.factorize(n)
    elif algorithm == 'rho_pollard':
        return factorization.rho_pollard.factorize(n, primality_tests=primality_tests)
    elif algorithm == 'p1_pollard':
        return factorization.p1_pollard.factorize(n, primality_tests=primality_tests)
    elif algorithm == 'garbell_quadratic':
        return factorization.garbell_quadratic.factorize(n, primality_tests=primality_tests)
    elif algorithm == 'garbell_cos_nombres':
        return factorization.garbell_cos_nombres.factorize(n, primality_tests=primality_tests)
    else:
        # if there was no algorithm math we raise a value error
        error_msg = 'Unrecognized factorization algorithm {}.'.format(algorithm)
        raise ValueError(error_msg)
