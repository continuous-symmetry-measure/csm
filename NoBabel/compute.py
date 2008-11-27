from math import floor

def factorial(n):
    res = 1;
    for i in xrange(1,n+1):
        res *= i;
    return res

def compute(n):
    temp = 0.0
    for i in xrange(0, floor(n/2) + 1):
        temp += 1.0/factorial(n - 2 * i)/2**i/factorial(i)
    return factorial(n) * temp
