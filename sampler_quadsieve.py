import numpy as np
from sympy import randprime
from timeit import default_timer as timer

from functools import reduce
from math import log, ceil
import copy
from typing import Union, List, Sequence, Tuple
import os
from numpy import array, where
from typing import Sequence, List, Tuple

def gf2(matrix: Sequence[Sequence[int]]) -> List[List[int]]:
    """Returns GF(2) form of matrix (only binary elements).

    :param matrix: matrix
    :return: matrix over GF(2)
    """
    return [[x % 2 for x in row] for row in matrix]

def fast_gauss(matrix: Sequence[Sequence[int]]) -> Tuple[List[List[int]], List[bool]]:
    """Performs a fast Gaussian Elimination over GF(2) on matrix.

    :param matrix: matrix to triangularize
    :return: triangularized matrix of zeros and ones.
    """
    m = array(gf2(matrix))

    marked = [False] * len(m)

    for j, column in enumerate(m.T):
        try:
            pivot = where(column == 1)[0][0]
            marked[pivot] = True

            for k, col in enumerate(m.T):
                if k == j:
                    continue

                if col[pivot] == 1:
                    m[:, k] += m[:, j]
                    m[:, k] %= 2

        except (ValueError, IndexError):
            pass

    return m.tolist(), marked

def perfect_square_combinations(base: Sequence[Sequence[int]]) -> List[List[int]]:
    """Generator for finding all perfect squares possible to form using factors from base.

    :param base: list of lists of exponents
    :return: list of lists containing indexes of base elements to use to form a perfect square
    """
    m, marked = fast_gauss(base)

    rows_independent = [r for i, r in enumerate(m) if marked[i]]
    rows_dependent = [(i, row) for i, row in enumerate(m) if not marked[i]]

    for i, row in rows_dependent:
        ones_cols = [j for j, x in enumerate(row) if x == 1]

        rows_reducing = []

        for c in ones_cols:
            # Get a whole column in independent rows and search for 1
            r = [r[c] for r in rows_independent].index(1)
            rows_reducing.append(r)

        yield rows_reducing + [i]


BASE_DIR = os.path.dirname(__file__)

PRIMES_URI = BASE_DIR + "/primes.txt"


def legendre(a: int, p: int) -> int:
    """Compute the Legendre symbol a|p using Euler's criterion.

    :param a:
    :param p: Odd prime.
    :return: 0 if a is divisible by b, 1 if a has a square root modulo p, 1 otherwise
    """

    if p % 2 == 0:
        raise ValueError('p should be an odd prime.')

    ls = pow(a, (p - 1) // 2, p)
    return -1 if ls == p - 1 else ls


def factorize(n: int, primes: Sequence[int]) -> Union[List[int], None]:
    """Factorizes n using only given primes.

    Returns exponent vector of n.
    If that's not possible, returns None.

    :param n: number to factorize
    :param primes: list of primes
    :return: exponent vector or None
    """
    exp = [0] * len(primes)

    if n == 0:
        return exp

    for i, p in enumerate(primes):
        while n % p == 0:
            exp[i] += 1
            n /= p

        if n <= 1:
            return exp

    return None

class Base:

    def __init__(self, n, primes, base_size):
        self.n = n
        self.width = len(primes)
        self.size = base_size
        self.primes = copy.deepcopy(primes)
        self.x = []
        self.x_sq_minus_n = []
        self.x_sq_minus_n_exp = []

        self.generate()

    def generate(self):
        x = ceil(self.n ** 0.5)

        i = 0

        while (i < self.size):
            t = x ** 2 - self.n
            t_exp = factorize(t, self.primes)

            if t_exp is not None:
                i += 1
                self.x.append(x)
                self.x_sq_minus_n.append(t)
                self.x_sq_minus_n_exp.append(t_exp)
            elif t > self.n:
                self.x.clear()
                break

            x += 1

def get_smooth_primes(b: int = None) -> List[int]:
    """Loads primes list from the built-in txt base.

    :param b: The smoothness bound of primes to load.
    :return: List of primes.
    """

    primes = []

    with open(PRIMES_URI, 'r') as file:
        for line in file:
            for s in line.split():
                if not s.isnumeric():
                    continue

                p = int(s)

                if b is not None and p > b:
                    return primes

                primes.append(int(s))

    return primes


def smoothness_bound(n: int, mode: int) -> int:
    """Computes the optimal smoothness bound for n.

    Based on "SMOOTH NUMBERS AND THE QUADRATIC SIEVE" by Carl Pomerance.
    If n is sufficiently small, set mode to 1 for a higher bound or 2 for the highest bound.

    :param n: number being factorized
    :param mode: how conservatively the bound is set
    :return: smoothness bound number
    """
    x = ceil(log(n) * log(log(n)))
    x = x ** 0.5
    if mode == 0: x = x*0.707 # sqrt(2)/2
    x = 2.71 ** x
    return ceil(x) + 1


def value(exp_vector: Sequence[int], primes: Sequence[int]) -> int:
    """Computes value of number represented by exponent vector over certain primes.

    :param exp_vector: list of exponents
    :param primes: primes list
    :return: value of represented number
    """
    if exp_vector == [0] * len(primes):
        return 0

    val = 1

    for prime, exp in zip(primes, exp_vector):
        val *= prime ** exp

    return val


def gcd(a: int, b: int) -> int:
    while b != 0:
        t = b
        b = a % b
        a = t

    return a


def congruent(a: int, b: int, n: int) -> bool:
    return a % n == b % n


def solve(base: Base):
    n = base.n

    for rows in perfect_square_combinations(base.x_sq_minus_n_exp):
        exp_vector = [0] * base.width

        for i in rows:
            exp_vector = [x + base.x_sq_minus_n_exp[i][j] for j, x in enumerate(exp_vector)]

        exp_vector = [x // 2 for x in exp_vector]

        v = value(exp_vector, base.primes)

        xes_prod = reduce(lambda x, y: x * y, [base.x[i] for i in rows])

        a = v % n
        b = xes_prod % n

        if not congruent(a, b, n) and not congruent(a, -b, n):
            factor = gcd(a - b, n)

            if factor != 1 and factor != n:
                return [factor, n // factor]

    print('Didn\'t find any nontrivial solution. Try changing the smoothness bound and base size.')
    return None


def qs_factorize(n: int, b: Union[int, None] = None, mode: int = 0, base_size: Union[int, None] = None,
              primes: Union[List[int], None] = None):
    """Finds factors of n using quadratic sieve algorithm with b-smooth base.

    :param n: number being factorized
    :param b: smoothness bound
    :param base_size: size of sieve's base
    :param primes: primes base to use
    :param loud_mode: True for displaying output information while work
    :return: tuple of n's factors or None
    """

    negative = False

    if n < 0:
        n *= -1
        negative = True

    if b is None:
        b = smoothness_bound(n, mode)

    print('Smoothness bound: ', b)

    if primes is None:
        print('Loading the primes base...')

        primes = get_smooth_primes(b)

        primes = [2] + [p for p in primes[1:] if legendre(n, p) == 1]

    print('Primes in the base: ', len(primes))

    if base_size is None:
        base_size = len(primes)

    print('Size of the base: ', base_size)
    print('Generating the base...')
    base = Base(n, primes, base_size)
    if len(base.x) > 0: print('The base has been generated.')
    else:
        print("Error generating base!")
        return None

    print('Searching for the right combination...')

    solution = solve(base)

    if negative:
        solution *= -1

    return solution


####################################################################################################################

powMin = 10
powMax = 60
casesPerPow = 5

# Sla resultaten per tweemacht op
runtimeByPow = []
# Sla elke resultaat op
resultMatrix = []
xAxis = range(powMin, powMax+1)

for pow2 in range(powMin, powMax+1):
    print("Now running test cases on n = "+str(pow2))
    pow2over4 = pow(2, np.floor(pow2/2)-1)
    pow2over2 = pow2over4*2
    parity = pow2 % 2

    runtimesThisPow = []
    casesThisPow = []
    for i in range(0, casesPerPow):
        # Genereer twee willekeurige priemgetallen (p,q) van correcte grootte
        p = randprime(pow2over4, pow2over2)
        q = randprime(pow2over4+pow2over4*parity, pow2over2+pow2over2*parity)
        print("Case "+str(i+1)+": p = "+str(p)+", q = "+str(q))
        N = p*q
        thisMode = 0
        if pow2 <= 30: thisMode = 1
        casesThisPow.append(N)

        start2 = timer()

        factors = qs_factorize(N, mode=thisMode)

        end2 = timer()
        if (factors == None):
            print("Error: Matrix was irreducible or base was not constructible!")
            runtimesThisPow.append(-1)
        elif (factors[0] == p or factors[1] == p):
            runtimesThisPow.append(end2-start2)
            print("Result: factored in "+str(end2-start2)+" seconds")
        else:
            print("Error: factors were incorrect!")
            runtimesThisPow.append(-1)
    
    # Sla gemiddelde looptijd op
    averageRuntime = 0
    successfulCases = 0
    for time in runtimesThisPow:
        if(time > 0):
            averageRuntime += time
            successfulCases += 1
    if (successfulCases > 0): averageRuntime /= successfulCases
    runtimeByPow.append(averageRuntime)

    # Sla specifieke resultaten op
    for i in range(0, len(casesThisPow)):
        caseTuple = (casesThisPow[i], runtimesThisPow[i])
        resultMatrix.append(caseTuple)

with open("resultsQS.txt", "w") as file:
    caseIndex = 0
    for power in range(powMin, powMax+1):
        ind = power-powMin
        file.write("Average for n = "+str(xAxis[ind])+": "+str(runtimeByPow[ind])+" seconds. Specific cases: \n")
        for i in range(1, casesPerPow+1):
            if(resultMatrix[caseIndex][1] <= 0): file.write("  Case "+str(i)+": N = "+str(resultMatrix[caseIndex][0])+", unable to factor. \n")
            else: file.write("  Case "+str(i)+": N = "+str(resultMatrix[caseIndex][0])+", factored in "+str(resultMatrix[caseIndex][1])+" seconds. \n")
            caseIndex += 1
        
    file.write("Average runtimes as list: "+str(runtimeByPow))