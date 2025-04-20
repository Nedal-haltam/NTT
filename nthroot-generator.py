import sympy

# n = 4
# p = sympy.prime(974) # 7681
n = 4
PrimeNumber = 5
possible = sympy.nthroot_mod(1, n, PrimeNumber, all_roots=True)

finals = []
print(possible)
for poss in possible:
    print(poss)
    foo = True
    for i in range(1, n, 1):
        if poss ** i % PrimeNumber == 1:
            print(poss ** i % PrimeNumber)
            foo = False
    if foo:
        finals.append(poss)

print(finals)