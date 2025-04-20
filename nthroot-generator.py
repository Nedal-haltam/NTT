import sympy

n = 4
index = 974
# n = 8
# index = 234511
PrimeNumber = sympy.prime(index)
while (PrimeNumber - 1) % n != 0:
    index += 1
    PrimeNumber = sympy.prime(index)
roots = sympy.nthroot_mod(1, n, PrimeNumber, all_roots=True)

finals = []
print(roots)
for root in roots:
    foo = True
    for i in range(1, n, 1):
        if root ** i % PrimeNumber == 1:
            foo = False
            break
    if foo:
        finals.append(root)

print(f'The prime number is : {PrimeNumber}')
print(f'the possible nth roots of unity are {finals}')


# sequence  
seq = range(1, 5) 
prime_no = PrimeNumber
# ntt 
transform = sympy.ntt(seq, prime_no) 
invtransform = sympy.intt(transform, prime_no)
print ("NTT : ", transform) 
print ("inverse NTT : ", invtransform) 