#include <iostream>
#include <algorithm>
#include <vector>
#include <complex>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <iomanip>

using namespace std;
const double PI = acos(-1);

const int mod = 998244353; // A commonly used NTT-friendly prime
const int root = 3;        // Primitive root modulo mod
int mod_pow(int base, int exp, int m) {
    int result = 1;
    while (exp > 0) {
        if (exp & 1) result = (1LL * result * base) % m;
        base = (1LL * base * base) % m;
        exp >>= 1;
    }
    return result;
}

void ntt(vector<int>& a, bool invert, bool other) {
    int n = a.size();
    int root_pw = mod_pow(root, (mod - 1) / n, mod);
    if (invert) root_pw = mod_pow(root_pw, mod - 2, mod); // inverse root

    // Bit-reversal permutation
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
            int lencopy = (int)pow(2, i);
            int wlen = mod_pow(root, (mod - 1) / lencopy, mod);
            if (invert) wlen = mod_pow(wlen, mod - 2, mod);
            for (int j = 0; j < n / lencopy; j++) {
                int w = 1;
                for (int k = 0; k < lencopy / 2; k++) {
                    int u = a[(j * lencopy) + k];
                    int v = (1LL * a[(j * lencopy)+ k + lencopy / 2] * w) % mod;
                    a[(j * lencopy)+ k] = (u + v) % mod;
                    a[(j * lencopy)+ k + lencopy / 2] = (u - v + mod) % mod;
                    w = (1LL * w * wlen) % mod;
                }
            }
        }
    }
    else
    {
        for (int len = 2; len <= n; len <<= 1) {
            int wlen = mod_pow(root, (mod - 1) / len, mod);
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

void fft(vector<complex<double>>& a, bool invert) {
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

    // Cooley�Tukey FFT
    for (int len = 2; len <= n; len <<= 1) {
        double angle = sign * 2 * PI / len;
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

    // Normalize if inverse
    if (invert) {
        for (complex<double>& x : a)
            x /= n;
    }
}

vector<complex<double>> FFTMultiply(const vector<complex<double>>& a, const vector<complex<double>>& b) {
    vector<complex<double>> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n = 1;
    while (n < a.size() + b.size())
        n <<= 1;
    fa.resize(n);
    fb.resize(n);

    fft(fa, false);
    fft(fb, false);
    for (int i = 0; i < n; ++i)
        fa[i] *= fb[i];
    fft(fa, true);

    return fa;
}

vector<int> NTTMultiply(vector<int> const& a, vector<int> const& b, bool other) {
    vector<int> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n = 1;
    while (n < a.size() + b.size()) n <<= 1;
    fa.resize(n);
    fb.resize(n);

    ntt(fa, false, other);
    ntt(fb, false, other);
    for (int i = 0; i < n; i++)
        fa[i] = (1LL * fa[i] * fb[i]) % mod;
    ntt(fa, true, other);

    return fa;
}

vector<int> Multiply(vector<int> const& a, vector<int> const& b)
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

void FillVectorRandom(vector<int>& vec, int a, int b)
{
    for(int i = 0; i < vec.size(); i++)
        vec[i] = a + (rand() % (b - a + 1));
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
int main()
{
    int n = 10 * 1024 * 1024;
    srand(time(0));

    vector<int> a(n);
    FillVectorRandom(a, 1, 5);
    vector<int> b(n);
    FillVectorRandom(b, 1, 5);

    vector<int> c;

    StartClock();
        c = NTTMultiply(a, b, false);
    EvaluateClock("the origianl");

    c.clear();

    StartClock();
        c = NTTMultiply(a, b, true);
    EvaluateClock("the new one with indexing");

    return 0;
}