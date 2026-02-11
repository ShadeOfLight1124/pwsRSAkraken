import numpy as np
from sympy import randprime
from timeit import default_timer as timer

def SieveOfAtkin(L: int):
    # Maak de lijst van markers. De positie van de marker is het getal waarvoor het bijhoudt of het priem is.
    isPrime = np.array(range(0, L+1))
    print("possible to make array")
    for i in range(0, L+1):
        isPrime[i] = False
    print("loop 1")
    
    # Maak 2 en 3 priem.
    isPrime[2] = True
    isPrime[3] = True

    x = 1
    # Bereken de kwadraten een keer om tijd te besparen
    x2 = 1
    # Controleer op criteria
    while x2 <= L:
        y = 1
        y2 = y*y
        while y2 <= L:
            # Bereken de kwadratische uitdrukkingen een keer om tijd te besparen.
            n1 = 4*x2+y2
            # Criterium 1: 4x^2+y^2 is equivalent aan 1 modulo 4 en niet aan 1 modulo 6.
            if(n1 <= 3):
                y2 += 2*y+1
                y += 1
                continue
            if(n1 <= L and (n1 % 12 == 1 or n1 % 12 == 5)):
                isPrime[n1] = not isPrime[n1]
            
            # Criterium 2: 3x^2+y^2 is equivalent aan 1 modulo 6 en niet aan 1 modulo 4.
            # Gebruik vooral additieve operaties om tijd te besparen. Vermenigvuldig zo veel mogelijk alleen met tweemachten, omdat dat efficient is.
            n2 = n1-x2
            if(n2 <= 3):
                y2 += 2*y+1
                y += 1
                continue
            if(n2 <= L and n2 % 12 == 7):
                isPrime[n2] = not isPrime[n2]
            
            # Criterium 3: 3x^2-y^2 is equivalent aan 11 modulo 12.
            n3 = n2-2*y2
            if(n3 <= 3 or x<=y):
                y2 += 2*y+1
                y += 1
                continue
            if(n3 <= L and n3 % 12 == 11):
                isPrime[n3] = not isPrime[n3]

            # (y+1)^2 = y^2+2y+1, dus deze snellere bewerkingen zijn equivalent aan eerst 1 bij y optellen en dan y^2 opnieuw berekenen.
            y2 += 2*y+1
            y += 1
        x2 += 2*x+1
        x += 1
    
    print("loop 2")

    # Markeer alle getallen deelbaar door kwadraten als composiet.
    i = 5
    i2 = i*i
    while i2 <= L:
        if isPrime[i]:
            for j in range(i*i, L+1, i*i):
                isPrime[j] = False
        i2 += 2*i+1
        i += 1
    
    primes = []
    for i in range(2, L+1):
        if(isPrime[i]): primes.append(i)
    print("loop 3")
    
    return primes



####################################################################################################################

powMin = 10
powMax = 56
casesPerPow = 5

# Sla resultaten per tweemacht op
runtimeByPow = []
sievetimeByPow = []
# Sla elke resultaat op
resultMatrix = []
xAxis = range(powMin, powMax+1)

for pow2 in range(powMin, powMax+1):
    print("Now running test cases on n = "+str(pow2))
    pow2over4 = int(pow(2, np.floor(pow2/2)-1))
    pow2over2 = int(pow2over4*2)
    parity = pow2 % 2

    start1 = timer()
    primeList = SieveOfAtkin(pow2over2+1)
    end1 = timer()
    sievetimeByPow.append(end1-start1)

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

        factors = []
        for prime in primeList:
            if (N % prime == 0):
                factors.append(prime)
                factors.append(N // prime)

        end2 = timer()
        if (len(factors) == 0):
            print("Error: no factors were found!")
            runtimesThisPow.append(-1)
        elif (factors[0] == p or factors[1] == p):
            runtimesThisPow.append(end2-start2)
            print("Result: factored in "+str(end2-start2)+" seconds")
    
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

with open("resultsAtkin.txt", "w") as file:
    caseIndex = 0
    for power in range(powMin, powMax+1):
        ind = power-powMin
        file.write("Average for n = "+str(xAxis[ind])+": "+str(runtimeByPow[ind])+" seconds. Sieve time: "+str(sievetimeByPow[ind])+". Specific cases: \n")
        for i in range(1, casesPerPow+1):
            if(resultMatrix[caseIndex][1] <= 0): file.write("  Case "+str(i)+": N = "+str(resultMatrix[caseIndex][0])+", unable to factor. \n")
            else: file.write("  Case "+str(i)+": N = "+str(resultMatrix[caseIndex][0])+", factored in "+str(resultMatrix[caseIndex][1])+" seconds. \n")
            caseIndex += 1
    
    file.write("Sieve times as list: "+str(sievetimeByPow)+"\n")
    file.write("Average runtimes as list: "+str(runtimeByPow))