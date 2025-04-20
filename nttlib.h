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

void ntt(vector<int>& a, bool invert, bool other, int root, int mod) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j)
        {
            swap(a[i], a[j]);
            //cout << i << " " << j << endl;
        }
    }
    if (other)
    {
        for (int i = 1; i <= (int)log2(n); i++) {
            int len = 1 << i;
            int wlen = root;
            if (invert) wlen = mod_pow(wlen, mod - 2, mod);
            for (int j = 0; j < n / len; j++) {
                int w = 1;
                for (int k = 0; k < len / 2; k++) {
                    int u = a[(j * len) + k];
                    int v = (1LL * a[(j * len)+ k + len / 2] * w) % mod;
                    a[(j * len)+ k] = (u + v) % mod;
                    a[(j * len)+ k + len / 2] = (u - v + mod) % mod;
                    w = (1LL * w * wlen) % mod;
                }
            }
        }
    }
    else
    {
        for (int len = 2; len <= n; len <<= 1) {
            int wlen = root;
            if (invert) wlen = mod_pow(wlen, mod - 2, mod);
            for (int i = 0; i < n; i += len) {
                int w = 1;
                for (int j = 0; j < len / 2; ++j) {
                    int u = a[i + j];
                    int v = (1LL * a[i + j + len / 2] * w) % mod;
                    a[i + j] = (u + v) % mod;
                    a[i + j + len / 2] = (u - v + mod) % mod;
                    w = (1LL * w * wlen) % mod;
                }
            }
        }
    }

    if (invert) {
        int n_inv = mod_pow(n, mod - 2, mod);
        for (int& x : a) x = (1LL * x * n_inv) % mod;
    }
}



void fft(vector<complex<double>>& a, bool invert, bool other) {
    int n = a.size();
    int sign = invert ? -1 : 1;
    // Bit-reversal permutation
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j)
        {
            swap(a[i], a[j]);
        }
    }

    if (other)
    {
        for (int i = 1; i <= (int)log2(n); i++) {
            int len = 1 << i;
            double angle = sign * 2 * M_PI / len;
            complex<double> wlen(cos(angle), sin(angle));
            for (int j = 0; j < n / len; j++) {
                complex<double> w(1);
                for (int k = 0; k < len / 2; k++) {
                    complex<double> u = a[(j * len) + k];
                    complex<double> v = a[(j * len) + k + len / 2] * w;
                    a[(j * len) + k] = u + v;
                    a[(j * len) + k + len / 2] = u - v;
                    w *= wlen;
                }
            }
        }
    }
    else
    {
        // Cooley�Tukey FFT
        for (int len = 2; len <= n; len <<= 1) {
            double angle = sign * 2 * M_PI / len;
            complex<double> wlen(cos(angle), sin(angle));
            for (int i = 0; i < n; i += len) {
                complex<double> w(1);
                for (int j = 0; j < len / 2; ++j) {
                    complex<double> u = a[i + j];
                    complex<double> v = a[i + j + len / 2] * w;
                    a[i + j] = u + v;
                    a[i + j + len / 2] = u - v;
                    w *= wlen;
                }
            }
        }
}

    // Normalize if inverse
    if (invert) {
        for (complex<double>& x : a)
            x /= n;
    }
}

vector<complex<double>> FFTMultiply(const vector<complex<double>>& a, const vector<complex<double>>& b, bool other) {
    vector<complex<double>> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n = 1;
    while (n < a.size() + b.size())
        n <<= 1;
    fa.resize(n);
    fb.resize(n);

    fft(fa, false, other);
    fft(fb, false, other);
    for (int i = 0; i < n; ++i)
        fa[i] *= fb[i];
    fft(fa, true, other);

    return fa;
}

vector<int> NTTMultiply(vector<int> const& a, vector<int> const& b, bool other, int root, int mod) {
    vector<int> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n = 1;
    while (n < a.size() + b.size()) n <<= 1;
    fa.resize(n);
    fb.resize(n);

    ntt(fa, false, other, root, mod);
    ntt(fb, false, other, root, mod);
    for (int i = 0; i < n; i++)
        fa[i] = (1LL * fa[i] * fb[i]) % mod;
    ntt(fa, true, other, root, mod);

    return fa;
}

vector<int> NaiveMultiply(vector<int> const& a, vector<int> const& b)
{
    int n = a.size();
    int m = b.size();
    vector<int> c(n + m - 1, 0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            c[i + j] += a[i] * b[j];
        }
    }
    return c;
}

template<class T>
void DisplayVector(vector<T> vec)
{
    for (T& e : vec)
        cout << e << " ";
    cout << endl;
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

void NTTMultiplyBenchMark(int n, int root, int mod)
{
    int lower(1);
    int upper(5);
    auto RadomFunction = GetRandomInt;
    vector<int> a(n);
    FillVectorRandom(a, lower, upper, RadomFunction);
    vector<int> b(n);
    FillVectorRandom(b, lower, upper, RadomFunction);
    
    StartClock();
    vector<int> c = NTTMultiply(a, b, true, root, mod);
    EvaluateClock("");
    
    DisplayVector(a);
    DisplayVector(b);
    DisplayVector(c);
}

void FFTMultiplyBenchMark(int n)
{
    complex<double>lower(1, 0);
    complex<double>upper(5, 0);
    auto RadomFunction = GetRandomComplex;
    vector<complex<double>> a(n);
    FillVectorRandom(a, lower, upper, RadomFunction);
    vector<complex<double>> b(n);
    FillVectorRandom(b, lower, upper, RadomFunction);

    StartClock();
    vector<complex<double>> c = FFTMultiply(a, b, true);
    EvaluateClock("");
    DisplayVector(a);
    DisplayVector(b);
    DisplayVector(c);
}