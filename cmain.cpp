#include <iostream>
#include <algorithm>
#include <vector>
#include <complex>
#include <math.h>
#include <assert.h>
#include <time.h>

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

void ntt(vector<int>& a, bool invert) {
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

vector<int> NTTMultiply(vector<int> const& a, vector<int> const& b) {
    vector<int> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n = 1;
    while (n < a.size() + b.size()) n <<= 1;
    fa.resize(n);
    fb.resize(n);

    ntt(fa, false);
    ntt(fb, false);
    for (int i = 0; i < n; i++)
        fa[i] = (1LL * fa[i] * fb[i]) % mod;
    ntt(fa, true);

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

void FillVectorRandom(vector<int> vec, int a, int b)
{
    for(int i = 0; i < vec.size(); i++)
        vec[i] = a + (rand() % (b - a + 1));
}


int n = 8;
void GoldenFoo()
{
    for (int i = 2; i <= n; i <<= 1) {
        cout << "index i = " << i << "\n";
        for (int j = 0; j < n; j += i) {
            cout << "\tindex j = " << j << "\n";
            for (int k = 0; k < i / 2; k++) {
                cout << "\t\tindex k = " << k << "\n";
            }
        }
    }
}

void foo()
{
    for (int i = 1; i <= (int)log2(n); i++) {
        cout << "index i = " << (int)pow(2, i) << "\n";
        int icopy = (int)pow(2, i);
        for (int j = 0; j < n / icopy; j++) {
            cout << "\tindex j = " << j * icopy << "\n";
            for (int k = 0; k < icopy / 2; k++) {
                cout << "\t\tindex k = " << k << "\n";
            }
        }
    }
}


int main()
{
    GoldenFoo();
    cout << "---------------------------------------------------------------------------------------------------------------------------------------\n";
    foo();
    return 0;
    /*
    1 2 4 1 3 
    1 4 3 2 3 
    1 6 15 25 26 29 23 9 9 0 0 0 0 0 0 0
    */
    int n = 5;
    srand(time(0));
    vector<int> a(n);
    a = {1, 2, 4, 1, 3};
    vector<int> b(n);
    b = {1, 4, 3, 2, 3};

    clock_t tStart = clock();
    vector<int> c = NTTMultiply(a, b);
    cout << "done\n";
    printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
    DisplayVector(a);
    DisplayVector(b);
    DisplayVector(c);
    return 0;
}