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
    for (size_t i = 0; i < vec.size() - 1; i++)
    {
        cout << vec[i] << ", ";
    }
    cout << vec[vec.size() - 1];
    cout << "]\n";
}

template<class T>
void FillVectorRandom(vector<T>& vec, T a, T b, T (*GetRandomT)(T, T))
{
    for(size_t i = 0; i < vec.size(); i++)
        vec[i] = GetRandomT(a, b);
}

template<class T>
bool CompareVectors(vector<T>& c, vector<T>& cother)
{
    assert(c.size() == cother.size());
    for (size_t i = 0; i < c.size(); i++)
    {
        if (c[i] != cother[i])
        {
            cout << "difference found\n";
            return false;
        }
    }
    return true;
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


void bit_reverse_permutation(std::vector<int>& data) {
    int n = data.size();
    int bits = std::log2(n);
    
    for (int i = 0; i < n; ++i) {
        int rev = 0;
        int x = i;
        for (int j = 0; j < bits; ++j) {
            rev = (rev << 1) | (x & 1);
            x >>= 1;
        }
        if (i < rev) {
            std::swap(data[i], data[rev]);
        }
    }
}

vector<int> ntt_naive(vector<int>& a, int pr, int mod, bool invert = false) {
    int n = a.size();
    vector<int> A(n, 0);
    int root = mod_pow(pr, (mod - 1) / n, mod);
    if (invert)
        root = mod_pow(root, mod - 2, mod);
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
        for (size_t i = 0; i < A.size(); i++)
        {
            A[i] = (1LL * A[i] * n_inv) % mod;    
        }
    }
    return A;
}

vector<int> NTTMultiplynaive(vector<int>& a, vector<int>& b, int root, int mod)
{
    vector<int> fa(a), fb(b);
    size_t n = 1;
    while (n < a.size() + b.size()) n <<= 1;
    fa.resize(n);
    fb.resize(n);
    fa = ntt_naive(fa, root, mod, false);
    fb = ntt_naive(fb, root, mod, false);
    for (size_t i = 0; i < n; i++)
    {
        fa[i] = (1LL * fa[i] * fb[i]) % mod;
    }
    return ntt_naive(fa, root, mod, true);
}

int main()
{
    srand(time(0));
    vector<int> a = vector<int>();
    int n = 16;
    int p = 7681;
    int pr = 17;
    for (int i = 0; i < n; i++)
        a.push_back(i+1);
    vector<int> ntt = ntt_naive(a, pr, p);
    vector<int> intt = ntt_naive(ntt, pr, p, true);
    cout << "The prime number is : " << p << "\n";
    cout << "The primitive root is : " << pr << "\n";
    DisplayVector(a);
    DisplayVector(ntt);
    DisplayVector(intt);

    return 0;
}



