import sympy.ntheory.residue_ntheory

n = 4
p = sympy.prime(974) # 7681
possible = sympy.nthroot_mod(1, n, p, all_roots=True)

finals = []
print(possible)
for poss in possible:
    print(poss)
    foo = True
    for i in range(1, n, 1):
        if poss ** i % p == 1:
            print(poss ** i % p)
            foo = False
    if foo:
        finals.append(poss)

print(finals)