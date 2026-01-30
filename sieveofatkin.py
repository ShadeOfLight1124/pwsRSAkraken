import numpy as np

def SieveOfAtkin(L):
    # Maak de lijst van markers. De positie van de marker is het getal waarvoor het bijhoudt of het priem is.
    isPrime = np.array(range(0, L+1))
    for i in range(0, L+1):
        isPrime[i] = False
    
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

    # Markeer alle getallen deelbaar door kwadraten als composiet.
    i = 5
    i2 = i*i
    while i2 <= L:
        if isPrime[i]:
            for j in range(i*i, L+1, i*i):
                isPrime[j] = False
        i2 += 2*i+1
        i += 1
    
    # Geef alle priemgetallen weer.
    for p in range(2, L+1):
        if isPrime[p]:
            print(p)

SieveOfAtkin(100000000)