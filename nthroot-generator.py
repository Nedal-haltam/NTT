import sympy

def GenerateNTT(n : int, index : int):
    PrimeNumber = sympy.prime(index)
    while not ((PrimeNumber - 1) % n == 0 and PrimeNumber % (2 * n) == 1):
        index += 1
        PrimeNumber = sympy.prime(index)
    roots = sympy.nthroot_mod(1, n, PrimeNumber, all_roots=True)

    finals = []
    # print(roots)
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
    seq = range(1, n + 1)
    prime_no = PrimeNumber
    # ntt 
    transform = sympy.ntt(seq, prime_no) 
    invtransform = sympy.intt(transform, prime_no)
    print("NTT :",transform) 
    # print ("inverse NTT :",invtransform) 
    print('----------------------------------------------------------------------------------------------------------------')


GenerateNTT(8, 235511)

