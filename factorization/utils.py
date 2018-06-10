from sage.all import *


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

        # check that the inputs are correct and store them.
        # value.
        self._is_valid_value(value)
        self._value = abs(value)
        # factors.
        if factors is None:
            self._factors = {}
        else:
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

        out_str = 'The integer {} can be decomposed as:\n'.format(number)

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

        out_str = out_str if self._reduced_value == 1 else '{}{}*'.format(out_str, self._reduced_value)

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
    if factorization1.reduced_value % factorization2.value != 0:
        raise ValueError('Second factorization value does not divide first factorization reduced value')

    # Extract a list of integers from the dictionary of factors.
    power_factors = factorization2.factors
    factors = []
    for factor in power_factors:
        factors.extend([factor] * power_factors[factor])

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


def get_power_mod_n(base, power, mod):
    """
    Applies an optimized algorithm to get base^power % mod in logarithmic time and with no overflow problems.
    :param base: base.
    :param power: power to be given to the base.
    :param mod: modulo.
    :return: base^power % mod.
    """
    # We apply an optimized algorithm to compute integer powers in logarithmic time.
    prod = 1
    ratio = base % mod  # do this to avoid overflow
    # while we write the power in binary as power=i_ni_{n-1}...i_1 where i_j are either 0 or 1
    # we compute the product as b^(i_0)(b^2)^(i_1)...(b^{n+1})^(i_n).
    while power != 0:
        if (power % 2) == 1:
            prod *= ratio
            prod = prod % mod  # do this to avoid overflow.
            power -= 1
        # increase the ration and divide the power by 2 to get the next binary digit.
        ratio *= ratio
        ratio = ratio % mod  # do this to avoid overflow.
        power /= 2

    return prod


def get_prod_mod_n(factor_list, mod):
    """
    :param factor_list: Factors to be multiplied.
    :param mod: modulo.
    :return: product of the integers in the factor list modulo mod.
    """
    # we multiply all elements in the factor list making modulo at every step to avoid numeric problems
    prod = 1
    for factor in factor_list:
        prod = (prod * factor) % mod  # do this to avoid overflow.

    return prod


def get_factorial_modulo(n, mod):
    """
    :param n: Number whose factorial we wil compute.
    :param mod: modulo.
    :return: n! % mod
    """
    prod = 1
    for i in range(1, n + 1):
        prod *= i
        prod = prod % mod
    return prod


def get_lower_primes(bound=1000):
    """
    Gets all primes lower than bound.
    To get these primes it applies the eratostenes sieve method.
    :param bound: a positive integer to be used as a bound. It must be greater than 1.
    :return: all primes lower than that bound.
    """
    # make sure the bound is positive and greater than 1
    bound = abs(bound)
    if bound <= 1:
        error_msg = 'Bound should be greater than 1 but it is {}.'.format(bound)
        raise ValueError(error_msg)

    # Initialize sieve.
    sieve = [i for i in range(2, bound + 1)]
    prime_list = []
    # apply_algorithm iteratively
    while len(sieve) != 0:
        # get first element of the sieve (it is a prime number).
        p = sieve.pop(0)
        prime_list.append(p)

        # remove all values of the sieve divisible by p.
        sieve = [i for i in sieve if i % p != 0]

    return prime_list


def least_absolute_residue(a, n):
    """
    Gets the least absolute residue of a modulo n (an integer between -n/2 and n/2 congruent to a modulo n).
    :param a: integer.
    :param n: modulo.
    :return: least absolute residue of a modulo n.
    """
    # Get residue.
    a = a % n
    # If it is lower or equal to n/2 we return it else we subtract n.
    return a if a <= n / 2 else (a - n)


def is_smooth(candidate, factor_base):
    """
    Checks if the candidate is smooth with respect to factor_base
    :param candidate: An integer that we want to check if is a good candidate or not.
    :param factor_base: A list of distinct primes and, possibly, -1.
    :return: True if the least absolute residue modulo mod of the square
    of candidate is smooth with respect to factor_base else False.
    """
    # Check if -1 is in the factor base and act accordingly.
    if -1 in factor_base:
        # remove -1 from factor base and take the absolute value of residue.
        factor_base = [factor for factor in factor_base if factor != -1]
        candidate = abs(candidate)
    else:
        # If residue is negative we return False since -1 in not in the factor_base.
        if candidate < 0:
            return False

    # Remove all factors on the factor base from the residue.
    for factor in factor_base:
        while candidate % factor == 0:
            candidate /= factor

    # If what remains is the number 1 we return True else False since there are more factors.
    return candidate == 1


def get_multiplicity(n, factor):
    """
    :param n: An integer.
    :param factor: A possible divisor of n. It must be different than 1 or 0.
    :return: A positive integer d such that factor^d divides n but factor^(d+1) doesn't.
    """
    # The case where factor == -1 we return 1 if n is negative and 0 else.
    if factor == -1:
        return int(n < 0)

    d = 0
    # get the multiplicity by dividing iteratively.
    while n % factor == 0:
        n /= factor
        d += 1

    return d


def get_vector_representation(n, factor_base):
    """
    :param n: An integer smooth with respect to the factor base.
    :param factor_base: A list of primes and, possibly, -1.
    :return: A vector representation of n in the factor base.
    """
    return [get_multiplicity(n, factor) for factor in factor_base]


def is_linear_dependent(vector_list):
    """
    :param vector_list: A list of lists of integers. The lists of integers will be treated as vectors.
    :return: A tuple consisting on a boolean and a list of zeros and ones.
             The boolean is True if, when writing the vectors as vectors in Z/(2) we obtain a linearly dependent set
             and false otherwise.
             If the boolean is true the list is a list of all zeros except for ones on a minimal linear dependent set.
             Else is None.
    """
    # Build the matrix in sagemath.
    height = len(vector_list)
    width = len(vector_list[0])
    m = matrix(Integers(2), height, width, vector_list)

    # Get the kernel.
    ker = m.left_kernel()
    # If the kernel is 0 we return the tuple (False, None),
    # else we return the tuple (True, v) where v is a non-zero element of the kernel.
    if ker.dimension() == 0:
        return False, None
    else:
        v = [int(i) for i in ker.basis()[0]]
        return True, v


def get_element_from_vector(vector, factor_base, mod):
    """
    :param vector: A list of integers corresponding to an integers expressed in a factor base.
    :param factor_base: A list of distinct primes and, possibly, -1.
    :param mod: A positive integer indicating a modulus. All operations are going to be performed modulus mod.
    :return: The integer corresponding to the vector modulus mod.
    """
    prod = 1
    for power, factor in zip(vector, factor_base):
        prod = (prod * get_power_mod_n(factor, power, mod)) % mod

    return prod


def get_order(a, n):
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
        prod = (prod * a) % n
        order += 1

    return 1


def get_jacobi_symbol(a, n):
    """
    :param a: a positive integer (we don't know how to compute fast the Jacobi symbol for negative integers).
    :param n: a positive odd integer.
    :return: the Jacobi symbol (a/n)
    """
    if a < 0:
        raise ValueError('a should be non negative.')
    if (n % 2) == 0:
        raise ValueError('n should be odd but it is not.')

    # get the residue of a modulo n.
    a = a % n
    # if a is zero then the Jacobi symbol is 0.
    if a == 0:
        return 0

    # Extract all powers of 2 multiplying a and compute the Jacobi symbol for those.
    prod = 1
    ratio = 1 if (n - 1) % 8 in [0, 6] else -1  # (2/n)=(-1)^((n^2-1)/8)
    while (a % 2) == 0:
        a /= 2
        prod *= ratio

    # if it is a power of 2 then the Jacobi symbol is prod.
    if a == 1:
        return prod

    # Now that a is odd we can apply the quadratic reciprocity.
    prod = prod if ((a - 1) % 4) == 0 or ((n - 1) % 4) == 0 else -prod  # we combine it with the previous results.
    return prod * get_jacobi_symbol(n, a)  # apply the algorithm recursively.
