import numpy as np
from sympy import randprime
from timeit import default_timer as timer


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
    pow2over4 = int(pow(2, np.floor(pow2/2)-1))
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
        casesThisPow.append(N)

        start2 = timer()

        factors = []
        for i in range(2, pow2over2+1):
            if (N % i == 0):
                factors.append(i)
                factors.append(N // i)
                break

        end2 = timer()

        if (factors[0] == p or factors[1] == p):
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

with open("resultsNaive.txt", "w") as file:
    caseIndex = 0
    for power in range(powMin, powMax+1):
        ind = power-powMin
        file.write("Average for n = "+str(xAxis[ind])+": "+str(runtimeByPow[ind])+" seconds. Specific cases: \n")
        for i in range(1, casesPerPow+1):
            if(resultMatrix[caseIndex][1] <= 0): file.write("  Case "+str(i)+": N = "+str(resultMatrix[caseIndex][0])+", unable to factor. \n")
            else: file.write("  Case "+str(i)+": N = "+str(resultMatrix[caseIndex][0])+", factored in "+str(resultMatrix[caseIndex][1])+" seconds. \n")
            caseIndex += 1
        
    file.write("Average runtimes as list: "+str(runtimeByPow))