# Elliptic Curve Division Polynomial Calculator

Interactive CLI that computes the **kernel polynomial** ψ_π(x) of an
endomorphism π = n + mφ of an elliptic curve

```
y² = x³ + Ax + B
```

for the two simplest cases of complex multiplication (CM):

| Case | Condition | Ring | Endomorphism φ | Action on points |
|---|---|---|---|---|
| Gaussian | `B = 0`, `A ≠ 0` | ℤ[i], i² = -1 | `[i]` | (x, y) ↦ (-x, iy) |
| Eisenstein | `A = 0`, `B ≠ 0` | ℤ[ω], ω = (-1+√-3)/2 | `[ω]` | (x, y) ↦ (ζ₃x, y) |

Given A, B, and π = n + mφ, the program builds ψ_π(x) out of the classical
division polynomials ψ_n, φ_n (Silverman, *The Arithmetic of Elliptic
Curves*, Ch. III.4) via the elliptic curve addition law, then factors and
LaTeX-prints the result. ψ_π(x) = 0 cuts out the x-coordinates of ker(π),
which is what makes this useful for isogeny/CM computations.

Two implementations are provided, doing the exact same computation with the
same recurrences — pick whichever runtime you have handy:

- **`dividingpolynomial.py`** — Python, using [SymPy](https://www.sympy.org/).
- **`dividingpolynomial.jl`** — Julia, using [SymPy.jl](https://github.com/JuliaPy/SymPy.jl)
  (a Julia wrapper around the same SymPy engine), so the two versions are
  verified to produce identical output.

`factor.py` is an unrelated scratch script (LaTeX → SymPy expression
parsing) and isn't part of this program.

## Usage

### Python

```bash
python3 -m venv venv
source venv/bin/activate
pip install sympy
python dividingpolynomial.py
```

### Julia

```bash
julia -e 'import Pkg; Pkg.add("SymPy")'
julia dividingpolynomial.jl
```

### Example session

```
Step 1: Define the Curve
   Enter A: 0
   Enter B: 1
   ✔ Detected CM by Eisenstein Integers ℤ[ω]

Step 2: Define the Map π = n + mω
   Enter n (Integer): 0
   Enter m (Complex): 3

KERNEL POLYNOMIAL Ψ_π(x)
(For π = 0 + 3ω)
------------------------------------------------------------
3x⁴ + 12x
------------------------------------------------------------
LATEX CODE:
\psi_{0+3ω}(x) = 3 x \left(x^{3} + 4\right)
```

## How it works

`get_division_polynomials(n, x, y, A, B)` memoizes the standard recurrence

```
ψ_0 = 0,  ψ_1 = 1,  ψ_2 = 2y
ψ_3 = 3x⁴ + 6Ax² + 12Bx - A²
ψ_4 = 4y(x⁶ + 5Ax⁴ + 20Bx³ - 5A²x² - 4ABx - 8B² - A³)

ψ_{2m}   = ψ_m (ψ_{m+2} ψ_{m-1}² - ψ_{m-2} ψ_{m+1}²) / (2y)
ψ_{2m+1} = ψ_{m+2} ψ_m³ - ψ_{m-1} ψ_{m+1}³

φ_0 = 1,  φ_1 = x,  φ_n = xψ_n² - ψ_{n+1}ψ_{n-1}
```

For a pure map (`n = 0` or `m = 0`) the kernel polynomial is just ψ_n or ψ_m.
For a mixed map `π = n + mφ`, it's built from ψ_n, φ_n, ψ_m, φ_m via the
group law, reduced modulo `y² = x³ + Ax + B`, then factored with
`sympy.factor`.
