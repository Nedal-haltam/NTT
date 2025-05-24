import sympy

def GenerateNTTFromPrime(n : int, PrimeNumber):
    print(f'The prime number is : {PrimeNumber}')
    print(f'The primitive root is : {sympy.primitive_root(PrimeNumber)}')
    seq = list(range(1, n + 1))
    # seq = [10, 913, 7679, 6764]
    transform = sympy.ntt(seq, PrimeNumber)
    print(seq)
    print(transform) 
    print (sympy.intt(transform, PrimeNumber)) 


n = 16
# p = sympy.prime(65421)
p = 7681
GenerateNTTFromPrime(n, p)

roots = sympy.nthroot_mod(1, n, p, all_roots=True)
validroots = []
for root in roots:
    a = pow(root, n, p) == 1
    b = True
    for i in range(1, n):
        if pow(root, i, p) == 1:
            b = False
            break
    if a and b:
        validroots.append(root)
print('roots: ', roots)
print('valid primitive n-th root of unity: ', validroots)
