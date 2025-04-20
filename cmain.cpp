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
            int power = (1LL * j * k) % n;
            int omega = mod_pow(root, power, mod);
            A[k] = (A[k] + a[j] * omega % mod) % mod;
        }
    }
    if (invert)
    {
        int n_inv = mod_pow(n, mod - 2, mod);
        for (int i = 0; i < A.size(); i++)
        {
            A[i] = (A[i] * n_inv) % mod;    
        }
    }
    return A;
}

int main()
{
    srand(time(0));

    int n = 4;
    int MOD = 7681;
    int ROOT = 3383;

    vector<int> a = {1, 2, 3, 4};
    // FillVectorRandom(a, 1, 5, GetRandomInt);
    vector<int> acopy(a);   

    DisplayVector(a);
    cout << "--------------------------------------------------------------------------------\n";
    ntt(a, false, false, ROOT, MOD);
    DisplayVector(a);
    cout << "--------------------------------------------------------------------------------\n";
    vector<int> c = ntt_naive(acopy, ROOT, MOD);
    DisplayVector(c);
    
    // NTTMultiplyBenchMark(n);
    // FFTMultiplyBenchMark(n);
    return 0;
}



