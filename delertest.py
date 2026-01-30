import numpy as np

# Delertest
def PrimalityByExhaustion(n): # returns bool
    d = 2
    # Bereken de kwadraat zonder vermenigvuldigingen.
    d2 = 4
    while (d2 <= n):
        if (n % d == 0): return d
        d2 += 2*d+1
        d += 1
    # n is een priem als het geen delers groter dan 1 en kleiner dan sqrt(n) heeft.
    return n

num = 7877*7919
factor = PrimalityByExhaustion(num)