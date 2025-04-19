package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand/v2"
	"time"
)

const PI = math.Pi

const mod int = 998244353 // A commonly used NTT-friendly prime
const root int = 3        // Primitive root modulo mod

func modPow(base, exp, mod int) int {
	result := 1
	base = base % mod
	for exp > 0 {
		if exp&1 == 1 {
			result = (result * base) % mod
		}
		base = (base * base) % mod
		exp >>= 1
	}
	return result
}

func ntt(a []int, invert bool) {
	n := len(a)
	rootPw := modPow(root, (mod-1)/n, mod)
	if invert {
		rootPw = modPow(rootPw, mod-2, mod) // inverse root
	}

	// Bit-reversal permutation
	for i, j := 1, 0; i < n; i++ {
		bit := n >> 1
		for ; j&bit != 0; bit >>= 1 {
			j ^= bit
		}
		j ^= bit
		if i < j {
			a[i], a[j] = a[j], a[i]
		}
	}

	for length := 2; length <= n; length <<= 1 {
		wlen := modPow(root, (mod-1)/length, mod)
		if invert {
			wlen = modPow(wlen, mod-2, mod)
		}
		for i := 0; i < n; i += length {
			w := 1
			for j := 0; j < length/2; j++ {
				u := a[i+j]
				v := (a[i+j+length/2] * w) % mod
				a[i+j] = (u + v) % mod
				a[i+j+length/2] = (u - v + mod) % mod
				w = (w * wlen) % mod
			}
		}
	}

	if invert {
		nInv := modPow(n, mod-2, mod)
		for i := range a {
			a[i] = (a[i] * nInv) % mod
		}
	}
}

func FFTMultiply(a, b []complex128) []complex128 {
	n := 1
	for n < len(a)+len(b) {
		n <<= 1
	}
	fa := make([]complex128, n)
	fb := make([]complex128, n)
	copy(fa, a)
	copy(fb, b)

	fft(fa, false)
	fft(fb, false)
	for i := 0; i < n; i++ {
		fa[i] *= fb[i]
	}
	fft(fa, true)

	return fa
}

func NTTMultiply(a, b []int) []int {
	n := 1
	for n < len(a)+len(b) {
		n <<= 1
	}
	fa := make([]int, n)
	fb := make([]int, n)
	copy(fa, a)
	copy(fb, b)

	ntt(fa, false)
	ntt(fb, false)
	for i := 0; i < n; i++ {
		fa[i] = (fa[i] * fb[i]) % mod
	}
	ntt(fa, true)

	return fa
}

func fft(a []complex128, invert bool) {
	n := len(a)
	sign := 1.0
	if invert {
		sign = -1.0
	}

	// Bit-reversal permutation
	for i, j := 1, 0; i < n; i++ {
		bit := n >> 1
		for ; j&bit != 0; bit >>= 1 {
			j ^= bit
		}
		j ^= bit
		if i < j {
			a[i], a[j] = a[j], a[i]
		}
	}

	// Cooley-Tukey FFT
	for length := 2; length <= n; length <<= 1 {
		angle := sign * 2 * PI / float64(length)
		wlen := cmplx.Rect(1, angle)
		for i := 0; i < n; i += length {
			w := complex(1, 0)
			for j := 0; j < length/2; j++ {
				u := a[i+j]
				v := a[i+j+length/2] * w
				a[i+j] = u + v
				a[i+j+length/2] = u - v
				w *= wlen
			}
		}
	}

	// Normalize if inverse
	if invert {
		for i := range a {
			a[i] /= complex(float64(n), 0)
		}
	}
}

func Multiply(a, b []int) []int {
	n := len(a)
	m := len(b)
	c := make([]int, n+m-1)

	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			c[i+j] += a[i] * b[j]
		}
	}

	return c
}

func main() {
	n := 2
	a := make([]int, n)
	b := make([]int, n)

	rand := rand.New(rand.NewPCG(uint64(time.Time.Nanosecond(time.Now())), uint64(time.Time.Nanosecond(time.Now()))))

	for i := 0; i < n; i++ {
		a[i] = 1 + rand.IntN(5)
		b[i] = 1 + rand.IntN(5)
	}

	start := time.Now()

	// Replace this with your FFT or NTT multiply function
	// c := FFTMultiply(a, b)
	fmt.Println("start")
	c := NTTMultiply(a, b)

	fmt.Println("done")
	fmt.Printf("Time taken: %.2fs\n", time.Since(start).Seconds())
	for i := range a {
		fmt.Print(a[i], " ")
	}
	fmt.Println()
	for i := range b {
		fmt.Print(b[i], " ")
	}
	fmt.Println()
	for i := range c {
		fmt.Print(c[i], " ")
	}
	fmt.Println()
}
