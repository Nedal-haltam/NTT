#include <iostream>
#include <algorithm>
#include <vector>
#include <complex>
#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>
#include <time.h>
#include <iomanip>

#include "nttlib.h"

using namespace std;







vector<int> ntt_naive(const vector<int>& a, int root, int mod, bool invert = false) {
    int n = a.size();
    vector<int> A(n, 0);
    int sign = (invert) ? -1 : 1;
    root = (invert) ? mod_pow(root, mod - 2, mod) : root;
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            int power = ((1LL * j * k) % n);
            int omega = mod_pow(root, power, mod);
            A[k] = (A[k] + 1LL * a[j] * omega % mod) % mod;
        }
    }
    if (invert)
    {
        int n_inv = mod_pow(n, mod - 2, mod);
        for (int i = 0; i < A.size(); i++)
        {
            A[i] = (1LL * A[i] * n_inv) % mod;    
        }
    }
    return A;
}

int main()
{
    srand(time(0));

    // int n = 4;
    // int PrimeMod = 7681;
    // int ROOT = 3383;
    int n = 8;
    int PrimeMod = 3263993;
    int ROOT = 1032801;

    // vector<int> a = {1, 2, 3, 4};
    vector<int> a = {1, 2, 3, 4, 5, 6, 7, 8};
    // FillVectorRandom(a, 1, 5, GetRandomInt);

    DisplayVector(a);
    vector<int> c = ntt_naive(a, ROOT, PrimeMod);
    DisplayVector(c);
    vector<int> c2 = ntt_naive(c, ROOT, PrimeMod, true);
    DisplayVector(c2);
    
    return 0;
}



