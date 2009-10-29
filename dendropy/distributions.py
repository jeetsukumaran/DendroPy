#! /usr/bin/env python

############################################################################
##  distribtions.py
##
##  Part of the DendroPy library for phylogenetic computing.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

"""
Methods random-variates and probabilities based on various parametric
distributions.
"""

import math
from dendropy import GLOBAL_RNG

def factorial(num):
    """factorial(n): return the factorial of the integer num.
    factorial(0) = 1
    factorial(n) with n<0 is -factorial(abs(n))
    """
    result = 1
    for i in xrange(1, abs(num)+1):
        result *= i
    return result

def binomial_coefficient(population, sample):
    "Returns  `population` choose `sample`."
    s = max(sample, population - sample)
    assert s <= population
    assert population > -1
    if s == population:
        return 1
    numerator = 1
    denominator = 1
    for i in xrange(s+1, population + 1):
        numerator *= i
        denominator *= (i - s)
    return numerator/denominator

def exp_pdf(value, rate):
    """
    Returns the probability density for an exponential distribution
    with an intensity of rate, evaluated at value.
    """
    return float(rate) * math.exp(-1.0 * rate * value)

def poisson_rv(rate, rng=None):
    """
    Returns a random number from a Poisson distribution with rate of
    `rate` (mean of 1/rate).
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
    process of `rate` over `period`.
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

def sample_multinomial(probs, rng=None):
    """Returns the index of the probability bin in `probs`.
    `probs` is assumed to sum to 1.0 (all rounding error contributes to the
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

def lengthed_choice(seq, lengths, rng=None):
    """
    Selects an element out of seq, with probabilities of each element
    given by the list `lengths` (which must be at least as long as the
    length of `seq` - 1).
    """
    if rng is None:
        rng = GLOBAL_RNG
    if lengths is None:
        lengths = [1.0/len(seq) for count in range(len(seq))]
    else:
        lengths = list(lengths)
    if len(lengths) < len(seq) - 1:
        raise Exception("Insufficient number of lengths specified")
    if len(lengths) == len(seq) - 1:
        lengths.append(1 - sum(lengths))
    prob_thresholds = []
    previous_break = 0.0
    for index in range(len(lengths)):
        prob_thresholds.append(previous_break + lengths[index])
        previous_break = prob_thresholds[index]
    pick = rng.random()
    for index, prob_threshold in enumerate(prob_thresholds):
        if pick <= prob_threshold:
            return seq[index]
    return seq[-1]
    
def chisqprob(chisq, df):
    """
    Returns the probability value associated with the provided chi-square
    value and df.  Adapted from chisq.c in Gary Perlman's |Stat.
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

def zprob(z):
    """
    Returns the probability value associated with the provided z-score.
    Adapted from z.c in Gary Perlman's |Stat.
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
    """ Geometric distribution per Devroye, Luc. _Non-Uniform Random Variate
    Generation_, 1986, p 500. http://cg.scs.carleton.ca/~luc/rnbookindex.html
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
    
    
