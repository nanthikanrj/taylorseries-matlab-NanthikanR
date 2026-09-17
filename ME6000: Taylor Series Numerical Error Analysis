# taylorseries-matlab-NanthikanR
ME 6000 Numerical Analysis  Assignment for Numerical Errors 2026-09-14
# Taylor Series Approximation for $e^x$
**Course:** ME 6000 Numerical Analysis  
**Assignment:** Numerical Errors & Floating-Point Arithmetic  
**Date:** September 14, 2026  
**Author:** Nanthikan Rutjirakul  

---

## 📌 Overview
This repository contains the MATLAB implementation for studying truncation and round-off errors in the Taylor series expansion of $e^x$:

$$e^x = \sum_{k=0}^{N} \frac{x^k}{k!} = 1 + x + \frac{x^2}{2!} + \frac{x^3}{3!} + \dots$$

The computations evaluate the effects of:
- **Precision:** Double-precision vs. Single-precision floating-point numbers.
- **Summation Order:** 
  - Increasing sequence: $k = 0, 1, 2, \dots, N$ (largest to smallest or vice versa depending on $x$)
  - Decreasing sequence: $k = N, N-1, \dots, 0$
- **Parameters:** Evaluated for $x = -10$ and $x = +10$ at orders $N = 10, 30, 50$.

---

## ⚙️ Implementation Details
To prevent arithmetic overflow during intermediate calculations (since $50! \approx 3.04 \times 10^{64}$ exceeds single-precision maximum capacity $\sim 3.40 \times 10^{38}$), terms are accumulated iteratively using the recurrence relation:

$$t_0 = 1, \quad t_k = t_{k-1} \cdot \frac{x}{k}$$

This guarantees all intermediate values remain well within the representable dynamic range while preserving true floating-point rounding behavior.

---

## 📂 Repository Contents
- **`taylor_exp_errors.m`**: Standard MATLAB script containing the full computation and helper functions.
- **`taylor_exp_errors.mlx`**: MATLAB Live Script format.
- **`taylor_exp_errors.pdf`**: Formatted document containing the code alongside complete execution tables.
