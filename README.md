# Number-Theoretic Transform (NTT) and Polynomial Multiplication

This project demonstrates and compares two approaches for polynomial multiplication:

* A **naive polynomial multiplication** implementation using direct coefficient-wise computation.
* A **Number-Theoretic Transform (NTT)** based multiplication for faster performance on large inputs over finite fields.

## Contents

### `cmain.cpp`

Implements the following:

* **Naive Polynomial Multiplication**
* **NTT-based Polynomial Multiplication** using a given prime modulus and primitive root
* Helper utilities for:

  * Vector display and comparison
  * Timing measurement
  * Random input generation
  * Bit-reversal permutation
  * Modular exponentiation

### `nthroot-generator.py`

* Uses `sympy` to:

  * Compute a primitive root modulo a given prime
  * Generate the NTT and its inverse on a sample sequence
  * Verify and display all valid `n`-th roots of unity modulo `p`

Ensure you have a C++ compiler that supports C++11 or later.

## Example Output of `cmain.cpp`

```
Naive NTT Polynomial Multiplication
The prime number is : 7681
The primitive root is : 17
The first polynomial is : 
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
The second polynomial is : 
[1, 2, 3, 4, 5]
The result of polynomial multiplication is : 
[1, 4, 10, 20, 35, 50, 65, 80, 95, 110, 114, 106, 85, 50, 0, 0]
--------------------------------
The prime number is : 7681
The primitive root is : 17
The original vector is : 
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
The NTT of the vector is : 
[136, 6734, 5637, 6199, 3652, 1868, 5998, 3125, 7673, 4540, 1667, 5797, 4013, 1466, 2028, 931]
The inverse NTT of the vector is : 
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
```

## Example Output of `nthroo-generator.py`

```
The prime number is : 7681
The primitive root is : 17
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
[136, 6734, 5637, 6199, 3652, 1868, 5998, 3125, 7673, 4540, 1667, 5797, 4013, 1466, 2028, 931]
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
roots:  [1, 527, 583, 849, 1213, 1728, 1925, 3383, 4298, 5756, 5953, 6468, 6832, 7098, 7154, 7680]
valid primitive n-th root of unity:  [527, 583, 849, 1728, 5953, 6832, 7098, 7154]
```

## Requirements

* **C++ Compiler** (e.g., `g++`)
* **Python 3** with `sympy` for root generation

## Notes

* `p = 7681` is a prime modulus suitable for NTT.
* `17` is a primitive root modulo `7681` and used for transform generation.
* The code includes correctness verification and runtime comparison capabilities.

---