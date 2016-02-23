#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Functions to calculate or draw values from various probability distributions.
"""

import math
from dendropy.calculate import combinatorics
from dendropy.utility import GLOBAL_RNG

def binomial_rv(n, p, rng=None):
    """
    Returns the number of successes in a sample of ``n`` trials, with the
    probability of success given by ``p``.  Using the BINV algorithm, as given
    by Kachitvicyanukul, V. and B. Schmeiser. 1988. Binomial random variate
    generation. Communications of the ACM 31: 216-222.
    Note: *NOT* the best algorithm according to the authors of the paper (who
    present their own as an alternative). Apart from rounding errors
    accumulating in the loop, it may also take a long time to return a value as
    ``n`` * ``p`` become large or even moderate (e.g., n=380 and p=0.8 still).
    """
    if rng is None:
        rng = GLOBAL_RNG
    q = 1 - p
    s = float(p) / q
    a = (n + 1) * s
    r = q ** n
    x = 0
    u = rng.random()
    while True:
        if u <= r:
            return x
        u = u - r
        x = x + 1
        r = (float(a)/x - s) * r

def exp_pdf(value, rate):
    """
    Returns the probability density for an exponential distribution
    with an intensity of rate, evaluated at value.
    """
    return float(rate) * math.exp(-1.0 * rate * value)

def poisson_rv(rate, rng=None):
    """
    Returns a random number from a Poisson distribution with rate of
    ``rate`` (mean of 1/rate).
    """
    if rng is None:
        rng = GLOBAL_RNG
    MAX_EXPECTATION = 64.0 # larger than this and we have underflow issues
    if rate > MAX_EXPECTATION:
        r = rate/2.0
        return poisson_rv(r) + poisson_rv(r)
    L = math.exp(-1.0 * rate)
    p = 1.0
    k = 0.0
    while p >= L:
        k = k + 1.0
        u = rng.random()
        p = p * u
    return int(k - 1.0)

def num_poisson_events(rate, period, rng=None):
    """
    Returns the number of events that have occurred in a Poisson
    process of ``rate`` over ``period``.
    """
    if rng is None:
        rng = GLOBAL_RNG
    events = 0
    while period > 0:
        time_to_next = rng.expovariate(1.0/rate)
        if time_to_next <= period:
            events = events + 1
        period = period - time_to_next
    return events

def poisson_pmf(k, rate):
    """
    Returns the probability of a number, ``k``, drawn from a Poisson distribution
    with rate parameter, ``rate`` (= 1/mean).
    """
    mean = 1.0/rate
    return float((mean ** k) * math.exp(-mean))/combinatorics.factorial(k)

def sample_multinomial(probs, rng=None):
    """Returns the index of the probability bin in ``probs``.
    ``probs`` is assumed to sum to 1.0 (all rounding error contributes to the
    last bin).
    """
    if rng is None:
        rng = GLOBAL_RNG
    u = rng.random()
    for n, i in enumerate(probs):
        u -= i
        if u < 0.0:
            return n
    return len(probs) - 1

def weighted_choice(seq, weights, rng=None):
    """
    Selects an element out of seq, with probabilities of each element
    given by the list ``weights`` (which must be at least as long as the
    length of ``seq`` - 1).
    """
    if rng is None:
        rng = GLOBAL_RNG
    if weights is None:
        weights = [1.0/len(seq) for count in range(len(seq))]
    else:
        weights = list(weights)
    if len(weights) < len(seq) - 1:
        raise Exception("Insufficient number of weights specified")
    if len(weights) == len(seq) - 1:
        weights.append(1 - sum(weights))
    return seq[weighted_index_choice(weights, rng)]

def weighted_index_choice(weights, rng=None):
    """
    (From: http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/)
    The following is a simple function to implement weighted random choice in
    Python. Given a list of weights, it returns an index randomly, according
    to these weights [1].
    For example, given [2, 3, 5] it returns 0 (the index of the first element)
    with probability 0.2, 1 with probability 0.3 and 2 with probability 0.5.
    The weights need not sum up to anything in particular, and can actually be
    arbitrary Python floating point numbers.
    If we manage to sort the weights in descending order before passing them
    to weighted_choice_sub, it will run even faster, since the random call
    returns a uniformly distributed value and larger chunks of the total
    weight will be skipped in the beginning.
    """
    if rng is None:
        rng = GLOBAL_RNG
    rnd = rng.random() * sum(weights)
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i

def chisq_pdf(chisq, df):
    """
    Returns the probability value associated with the provided chi-square
    value and df.  Adapted from chisq.c in Gary Perlman's Stat.
    """

    BIG = 20.0
    def ex(x):
        BIG = 20.0
        if x < -BIG:
            return 0.0
        else:
            return math.exp(x)

    if chisq <=0 or df < 1:
        return 1.0
    a = 0.5 * chisq
    if df%2 == 0:
        even = 1
    else:
        even = 0
    if df > 1:
        y = ex(-a)
    if even:
        s = y
    else:
        s = 2.0 * zprob(-math.sqrt(chisq))
    if (df > 2):
        chisq = 0.5 * (df - 1.0)
        if even:
            z = 1.0
        else:
            z = 0.5
        if a > BIG:
            if even:
                e = 0.0
            else:
                e = math.log(math.sqrt(math.pi))
            c = math.log(a)
            while (z <= chisq):
                e = math.log(z) + e
                s = s + ex(c*z-a-e)
                z = z + 1.0
            return s
        else:
            if even:
                e = 1.0
            else:
                e = 1.0 / math.sqrt(math.pi) / math.sqrt(a)
            c = 0.0
            while (z <= chisq):
                e = e * (a/float(z))
                c = c + e
                z = z + 1.0
            return (c*y+s)
    else:
        return s

def z_pmf(z):
    """
    Returns the probability value associated with the provided z-score.
    Adapted from z.c in Gary Perlman's Stat.
    """

    Z_MAX = 6.0    # maximum meaningful z-value
    if z == 0.0:
        x = 0.0
    else:
        y = 0.5 * math.fabs(z)
        if y >= (Z_MAX*0.5):
            x = 1.0
        elif (y < 1.0):
            w = y*y
            x = ((((((((0.000124818987 * w
                        -0.001075204047) * w +0.005198775019) * w
                      -0.019198292004) * w +0.059054035642) * w
                    -0.151968751364) * w +0.319152932694) * w
                  -0.531923007300) * w +0.797884560593) * y * 2.0
        else:
            y = y - 2.0
            x = (((((((((((((-0.000045255659 * y
                             +0.000152529290) * y -0.000019538132) * y
                           -0.000676904986) * y +0.001390604284) * y
                         -0.000794620820) * y -0.002034254874) * y
                       +0.006549791214) * y -0.010557625006) * y
                     +0.011630447319) * y -0.009279453341) * y
                   +0.005353579108) * y -0.002141268741) * y
                 +0.000535310849) * y +0.999936657524
    if z > 0.0:
                prob = ((x+1.0)*0.5)
    else:
                prob = ((1.0-x)*0.5)
    return prob


def geometric_rv(p, rng=None):
    """Geometric distribution per Devroye, Luc. Non-Uniform Random Variate
    Generation, 1986, p 500. http://cg.scs.carleton.ca/~luc/rnbookindex.html
    """
    if rng is None:
        rng = GLOBAL_RNG
    # p should be in (0.0, 1.0].
    if p <= 0.0 or p > 1.0:
        raise ValueError("p = %s: p must be in the interval (0.0, 1.0]" % p)
    elif p == 1.0:
        # If p is exactly 1.0, then the only possible generated value is 1.
        # Recognizing this case early means that we can avoid a log(0.0) later.
        # The exact floating point comparison should be fine. log(eps) works just
        # dandy.
        return 1

    # random() returns a number in [0, 1). The log() function does not
    # like 0.
    U = 1.0 - rng.random()

    # Find the corresponding geometric variate by inverting the uniform variate.
    G = int(math.ceil(math.log(U) / math.log(1.0 - p)))
    return G

def hypergeometric_pmf(x, m, n, k):
    """
    Given a population consisting of ``m`` items of class M and ``n`` items of class N,
    this returns the probability of observing ``x`` items of class M when sampling
    ``k`` times without replacement from the entire population (i.e., {M,N})

            p(x) = (choose(m, x) * choose(n, k-x)) / choose(m+n, k)
    """
    return float(combinatorics.choose(m, x) * combinatorics.choose(n, k-x))/combinatorics.choose(m+n, k)

def hypergeometric_pmf(x, m, n, k):
    """
    Given a population consisting of ``m`` items of class M and ``n`` items of class N,
    this returns the probability of observing ``x`` items of class M when sampling
    ``k`` times without replacement from the entire population (i.e., {M,N})

            p(x) = (choose(m, x) * choose(n, k-x)) / choose(m+n, k)
    """
    # following fails with 'OverflowError: long int too large to convert to
    # float' with large numbers
    # return float(combinatorics.choose(m, x) * combinatorics.choose(n, k-x))/combinatorics.choose(m+n, k)
    a = math.log(combinatorics.choose(m, x))
    b = math.log(combinatorics.choose(n, k-x))
    c = math.log(combinatorics.choose(m+n, k))
    return math.exp(a+b-c)
