from typing import List, Callable, Optional

from factorization.utils import Factorization
from factorization.primality_tests import trial_division


def factorize(n, primality_tests: Optional[List[Callable[[int], bool]]]=None):
    if primality_tests is None:
        primality_tests = [trial_division]
    factorization = Factorization(n)
    return factorization
