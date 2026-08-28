"""
Elliptic Curve Division Polynomial Calculator
==============================================

Computes the "kernel polynomial" psi_pi(x) that cuts out the kernel of an
endomorphism pi = n + m*s of an elliptic curve y^2 = x^3 + Ax + B that has
complex multiplication (CM) by either:

  * the Gaussian integers  Z[i]      (A != 0, B == 0), where s = i acts by
    (x, y) |-> (-x, i*y), or
  * the Eisenstein integers Z[omega] (A == 0, B != 0), where s = omega acts
    by (x, y) |-> (zeta_3 * x, y), zeta_3 a primitive cube root of unity.

The classical division polynomials psi_n and phi_n (Silverman, "The
Arithmetic of Elliptic Curves", Ch. III.4) satisfy the recurrences below and
determine the multiplication-by-n map. For a CM map pi = n + m*s, the kernel
polynomial is built out of psi_n, phi_n, psi_m, phi_m using the addition law
for pi = [n] + [m*s] on the curve.
"""

import sympy as sp


# --- TERMINAL STYLING CONSTANTS ---
class Style:
    """ANSI escape codes used to colorize the CLI output."""
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    END = '\033[0m'


def to_superscript(n):
    """Convert a non-negative integer to its unicode superscript form (e.g. 12 -> '¹²')."""
    digits = "0123456789"
    supers = "⁰¹²³⁴⁵⁶⁷⁸⁹"
    return "".join(supers[int(d)] for d in str(n))


def format_poly_clean(poly, var_char='x'):
    """
    Render a single-variable SymPy polynomial as one readable line, e.g.
    'x³ - 4x + 1', instead of SymPy's default '+ -4*x' style output.
    """
    poly = sp.Poly(poly, sp.symbols(var_char))
    term_list = poly.terms()

    if not term_list:
        return "0"

    terms = []
    for power_tuple, coeff in term_list:
        power = power_tuple[0]

        # Leading term has no '+', later terms get ' + ' / ' - ' depending on sign.
        if not terms:
            sign = "-" if coeff < 0 else ""
        else:
            sign = " - " if coeff < 0 else " + "

        abs_coeff = abs(coeff)

        # Omit an explicit coefficient of 1 unless the term is the constant.
        if abs_coeff == 1 and power != 0:
            coeff_str = ""
        else:
            coeff_str = str(abs_coeff)

        # x^0 -> "", x^1 -> "x", x^n -> "x" + superscript(n).
        if power == 0:
            var_str = ""
        elif power == 1:
            var_str = var_char
        else:
            var_str = f"{var_char}{to_superscript(power)}"

        terms.append(f"{sign}{coeff_str}{var_str}")

    return "".join(terms)


def get_division_polynomials(n_val, x, y, A, B):
    """
    Return (psi_n, phi_n), the n-th division polynomial and the associated
    x-coordinate numerator polynomial, for the curve y^2 = x^3 + Ax + B.

    psi_n satisfies the standard recurrence (Silverman III.4):
        psi_0 = 0, psi_1 = 1, psi_2 = 2y,
        psi_3 = 3x^4 + 6Ax^2 + 12Bx - A^2,
        psi_4 = 4y(x^6 + 5Ax^4 + 20Bx^3 - 5A^2x^2 - 4ABx - 8B^2 - A^3),
    and for m >= 2:
        psi_{2m}   = psi_m * (psi_{m+2} psi_{m-1}^2 - psi_{m-2} psi_{m+1}^2) / (2y)
        psi_{2m+1} = psi_{m+2} psi_m^3 - psi_{m-1} psi_{m+1}^3
    phi_n is the numerator of the x-coordinate of the n-torsion multiple:
        phi_0 = 1, phi_1 = x, phi_n = x*psi_n^2 - psi_{n+1}*psi_{n-1}.

    Results are memoized per call since the recurrence revisits the same
    index many times (e.g. computing psi_n typically needs psi_{n-2..n+2}).
    """
    psi_cache = {}

    def psi(n):
        if n in psi_cache:
            return psi_cache[n]
        if n == 0:
            return 0
        if n == 1:
            return 1
        if n == 2:
            return 2 * y
        if n == 3:
            return 3 * x**4 + 6 * A * x**2 + 12 * B * x - A**2
        if n == 4:
            return 4 * y * (x**6 + 5 * A * x**4 + 20 * B * x**3
                             - 5 * A**2 * x**2 - 4 * A * B * x - 8 * B**2 - A**3)

        m = n // 2
        if n % 2 == 0:
            res = (psi(m) * (psi(m + 2) * psi(m - 1)**2 - psi(m - 2) * psi(m + 1)**2)) / (2 * y)
        else:
            res = psi(m + 2) * psi(m)**3 - psi(m - 1) * psi(m + 1)**3
        psi_cache[n] = sp.expand(res)
        return res

    def phi(n):
        if n == 0:
            return 1
        if n == 1:
            return x
        return sp.expand(x * psi(n)**2 - psi(n + 1) * psi(n - 1))

    return psi(n_val), phi(n_val)


def main():
    print("\n" * 2)
    print(f"{Style.HEADER}{Style.BOLD}╔════════════════════════════════════════════════════╗")
    print(f"║   ELLIPTIC CURVE DIVISION POLYNOMIAL CALCULATOR    ║")
    print(f"║           For Complex Multiplication               ║")
    print(f"╚════════════════════════════════════════════════════╝{Style.END}")

    # --- STEP 1: DEFINE CURVE ---
    print(f"\n{Style.BOLD}Step 1: Define the Curve{Style.END}")
    print(f"Format: y² = x³ + {Style.CYAN}A{Style.END}x + {Style.CYAN}B{Style.END}")

    try:
        A_in = int(input(f"   Enter A: {Style.CYAN}"))
        print(f"{Style.END}", end="")
        B_in = int(input(f"   Enter B: {Style.CYAN}"))
        print(f"{Style.END}", end="")
    except ValueError:
        return

    # A curve is singular (not elliptic) iff its discriminant vanishes.
    disc = 4 * A_in**3 + 27 * B_in**2
    if disc == 0:
        print(f"\n{Style.RED}[CRITICAL ERROR] Singular Curve.{Style.END}")
        return

    cm_type = None
    s_symbol = ""

    # --- CM DETECTION & DESCRIPTION BLOCK ---
    # Only the two simplest CM cases are handled: j = 1728 (Z[i], B = 0)
    # and j = 0 (Z[omega], A = 0).
    if B_in == 0 and A_in != 0:
        cm_type = "Gaussian"
        s_symbol = "i"
        print(f"\n   {Style.GREEN}✔ Detected CM by Gaussian Integers ℤ[i]{Style.END}")
        print(f"   {Style.DIM}----------------------------------------{Style.END}")
        print(f"   • {Style.BOLD}Ring:{Style.END}   ℤ[i] where i² = -1")
        print(f"   • {Style.BOLD}Map s:{Style.END}  Corresponds to [i]")
        print(f"   • {Style.BOLD}Action:{Style.END} (x, y) ↦ (-x, iy)")
        print(f"   {Style.DIM}----------------------------------------{Style.END}")

    elif A_in == 0 and B_in != 0:
        cm_type = "Eisenstein"
        s_symbol = "ω"
        print(f"\n   {Style.GREEN}✔ Detected CM by Eisenstein Integers ℤ[ω]{Style.END}")
        print(f"   {Style.DIM}----------------------------------------{Style.END}")
        print(f"   • {Style.BOLD}Ring:{Style.END}   ℤ[ω] where ω = (-1 + √-3)/2")
        print(f"   • {Style.BOLD}Map s:{Style.END}  Corresponds to [ω]")
        print(f"   • {Style.BOLD}Action:{Style.END} (x, y) ↦ (ζ₃x, y)")
        print(f"           (where ζ₃ is a primitive cube root of unity)")
        print(f"   {Style.DIM}----------------------------------------{Style.END}")

    else:
        print(f"\n{Style.YELLOW}⚠ Warning: No simple CM detected.{Style.END}")
        return

    # --- STEP 2: DEFINE MAP ---
    print(f"\n{Style.BOLD}Step 2: Define the Map π = n + m{s_symbol}{Style.END}")
    try:
        n_map = int(input(f"   Enter n (Integer): {Style.CYAN}"))
        print(f"{Style.END}", end="")
        m_map = int(input(f"   Enter m (Complex): {Style.CYAN}"))
        print(f"{Style.END}", end="")
    except ValueError:
        return

    # Computation
    x, y = sp.symbols('x y')
    curve_eq = x**3 + A_in * x + B_in

    def clean_poly(expr):
        """Simplify and reduce modulo the curve equation y^2 = x^3 + Ax + B."""
        return sp.simplify(expr).subs(y**2, curve_eq)

    print(f"\n{Style.BLUE} ➤ Computing...{Style.END}")

    # --- EDGE CASE HANDLING ---
    if m_map == 0:
        # Case 1: Pure Integer Map [n] -- kernel is just psi_n.
        final_poly, _ = get_division_polynomials(abs(n_map), x, y, A_in, B_in)
        final_eq = clean_poly(final_poly)

    elif n_map == 0:
        # Case 2: Pure Complex Map [m]s -- kernel is psi_m.
        final_poly, _ = get_division_polynomials(abs(m_map), x, y, A_in, B_in)
        final_eq = clean_poly(final_poly)

    else:
        # Case 3: Mixed Map [n] + [m]s -- combine via the addition law.
        psi_n, phi_n = get_division_polynomials(abs(n_map), x, y, A_in, B_in)
        psi_m, phi_m = get_division_polynomials(abs(m_map), x, y, A_in, B_in)

        if cm_type == "Gaussian":
            term1 = clean_poly(phi_n * psi_m**2)
            term2 = clean_poly(phi_m * psi_n**2)
            final_eq = term1 + term2
        elif cm_type == "Eisenstein":
            zeta = sp.exp(2 * sp.pi * sp.I / 3)
            term1 = clean_poly(phi_n * psi_m**2)
            term2 = clean_poly(phi_m * psi_n**2)
            final_eq = term1 - zeta * term2

    result = sp.factor(final_eq)

    # --- OUTPUT ---
    print("\n" + "=" * 60)
    print(f"{Style.BOLD}KERNEL POLYNOMIAL Ψ_π(x){Style.END}")
    print(f"(For π = {n_map} + {m_map}{s_symbol})")
    print("-" * 60)

    try:
        clean_str = format_poly_clean(result)
        print(f"\n{Style.GREEN}{Style.BOLD}{clean_str}{Style.END}\n")
    except Exception:
        # format_poly_clean assumes a single-variable polynomial in x;
        # fall back to SymPy's own pretty-printer if that assumption fails.
        sp.init_printing(use_unicode=True)
        sp.pprint(result)

    print("-" * 60)
    print(f"{Style.YELLOW}LATEX CODE:{Style.END}")
    print(f"\\psi_{{{n_map}+{m_map}{s_symbol}}}(x) = " + sp.latex(result))
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
