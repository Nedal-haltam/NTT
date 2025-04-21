#include <iostream>
#include <algorithm>
#include <vector>
#include <complex>
#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>
#include <time.h>
#include <iomanip>

using namespace std;


int mod_pow(int base, int exp, int mod) {
    int result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1)
            result = (1LL * result * base) % mod;
        base = (1LL * base * base) % mod;
        exp >>= 1;
    }
    return result;
}

template<class T>
void DisplayVector(vector<T> vec)
{
    cout << "[";
    for (int i = 0; i < vec.size() - 1; i++)
    {
        cout << vec[i] << ", ";
    }
    cout << vec[vec.size() - 1];
    cout << "]\n";
}

template<class T>
void FillVectorRandom(vector<T>& vec, T a, T b, T (*GetRandomT)(T, T))
{
    for(int i = 0; i < vec.size(); i++)
        vec[i] = GetRandomT(a, b);
}

template<class T>
void CompareVectors(vector<T>& c, vector<T>& cother)
{
    assert(c.size() == cother.size());
    for (int i = 0; i < c.size(); i++)
    {
        if (c[i] != cother[i])
        {
            cout << "difference found\n";
            exit(EXIT_FAILURE);
        }
    }
}


clock_t clocker;
void StartClock()
{
    clocker = clock();
}
void EvaluateClock(string Message)
{
    clock_t t = clock() - clocker;
    cout << Message << "\n";
    cout << std::fixed << std::setprecision(3) << "Time taken: " << (double)(t) / CLOCKS_PER_SEC << "s\n";
}

int GetRandomInt(int a, int b)
{
    return a + (rand() % (b - a + 1));
}

complex<double> GetRandomComplex(complex<double> a, complex<double> b)
{
    complex<double> t(b - a + complex<double>(1, 1));
    return complex<double>(a.real() + rand() % (int)t.real(), a.imag() + rand() % (int)t.imag());
}


vector<int> ntt_naive(vector<int>& a, int root, int mod, bool invert = false) {
    int n = a.size();
    vector<int> A(n, 0);
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

vector<int> NTTMultiplynaive(vector<int>& a, vector<int>& b, int root, int mod)
{
    vector<int> fa(a), fb(b);
    int n = 1;
    while (n < a.size() + b.size()) n <<= 1;
    fa.resize(n);
    fb.resize(n);
    fa = ntt_naive(fa, root, mod, false);
    fb = ntt_naive(fb, root, mod, false);
    for (int i = 0; i < n; i++)
    {
        fa[i] = (1LL * fa[i] * fb[i]) % mod;
    }
    return ntt_naive(fa, root, mod, true);
}

int main()
{
    srand(time(0));
    vector<int> p1 = {1, 2, 3, 4};
    vector<int> p2 = {1, 2, 3, 4};

    int mod = 3278753;
    int root = 449739;

    vector<int> c = NTTMultiplynaive(p1, p2, root, mod);
    DisplayVector(c);

    return 0;
}



