import sympy as sp
import sys

# --- TERMINAL STYLING CONSTANTS ---
class Style:
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
    """Converts an integer to a unicode superscript string."""
    digits = "0123456789"
    supers = "⁰¹²³⁴⁵⁶⁷⁸⁹"
    return "".join(supers[int(d)] for d in str(n))

def format_poly_clean(poly, var_char='x'):
    """Manually formats a SymPy polynomial to a single clean line."""
    poly = sp.Poly(poly, sp.symbols(var_char))
    term_list = poly.terms()
    
    if not term_list: return "0"
    
    terms = []
    for power_tuple, coeff in term_list:
        power = power_tuple[0]
        
        # Handle Sign
        if not terms: sign = "-" if coeff < 0 else ""
        else: sign = " - " if coeff < 0 else " + "
        
        abs_coeff = abs(coeff)
        
        # Handle Coefficient
        if abs_coeff == 1 and power != 0: coeff_str = ""
        else: coeff_str = str(abs_coeff)
            
        # Handle Variable
        if power == 0: var_str = ""
        elif power == 1: var_str = var_char
        else: var_str = f"{var_char}{to_superscript(power)}"
            
        terms.append(f"{sign}{coeff_str}{var_str}")
        
    return "".join(terms)

def get_division_polynomials(n_val, x, y, A, B):
    psi_cache = {}
    
    def psi(n):
        if n in psi_cache: return psi_cache[n]
        if n == 0: return 0
        if n == 1: return 1
        if n == 2: return 2*y
        if n == 3: return 3*x**4 + 6*A*x**2 + 12*B*x - A**2
        if n == 4: return 4*y*(x**6 + 5*A*x**4 + 20*B*x**3 - 5*A**2*x**2 - 4*A*B*x - 8*B**2 - A**3)
        
        m = n // 2
        if n % 2 == 0:
            res = (psi(m) * (psi(m+2)*psi(m-1)**2 - psi(m-2)*psi(m+1)**2)) / (2*y)
        else:
            res = psi(m+2)*psi(m)**3 - psi(m-1)*psi(m+1)**3
        psi_cache[n] = sp.expand(res)
        return res

    def phi(n):
        if n == 0: return 1
        if n == 1: return x
        return sp.expand(x * psi(n)**2 - psi(n+1) * psi(n-1))

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

    # CM Check
    disc = 4*A_in**3 + 27*B_in**2
    if disc == 0:
        print(f"\n{Style.RED}[CRITICAL ERROR] Singular Curve.{Style.END}")
        return

    cm_type = None
    s_symbol = ""
    
    # --- CM DETECTION & DESCRIPTION BLOCK ---
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
    curve_eq = x**3 + A_in*x + B_in
    
    def clean_poly(expr):
        return sp.simplify(expr).subs(y**2, curve_eq)

    print(f"\n{Style.BLUE} ➤ Computing...{Style.END}")

    # --- EDGE CASE HANDLING ---
    if m_map == 0:
        # Case 1: Pure Integer Map [n]
        final_poly, _ = get_division_polynomials(abs(n_map), x, y, A_in, B_in)
        final_eq = clean_poly(final_poly)

    elif n_map == 0:
        # Case 2: Pure Complex Map [m]s
        final_poly, _ = get_division_polynomials(abs(m_map), x, y, A_in, B_in)
        final_eq = clean_poly(final_poly)

    else:
        # Case 3: Mixed Map [n] + [m]s
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
    print("\n" + "="*60)
    print(f"{Style.BOLD}KERNEL POLYNOMIAL Ψ_π(x){Style.END}")
    print(f"(For π = {n_map} + {m_map}{s_symbol})")
    print("-"*(60))
    
    try:
        clean_str = format_poly_clean(result)
        print(f"\n{Style.GREEN}{Style.BOLD}{clean_str}{Style.END}\n")
    except:
        sp.init_printing(use_unicode=True)
        sp.pprint(result)

    print("-"*(60))
    print(f"{Style.YELLOW}LATEX CODE:{Style.END}")
    print(f"\\psi_{{{n_map}+{m_map}{s_symbol}}}(x) = " + sp.latex(result))
    print("="*60 + "\n")

if __name__ == "__main__":
    main()
