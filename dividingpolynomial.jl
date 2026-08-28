#=
Elliptic Curve Division Polynomial Calculator (Julia port)
============================================================

Same program as dividingpolynomial.py, translated to Julia using SymPy.jl
(a Julia wrapper around Python's SymPy) so the two versions share identical
symbolic computations. See dividingpolynomial.py's module docstring for the
math background on division polynomials and the CM kernel-polynomial
construction.

Run with:  julia dividingpolynomial.jl
=#

using SymPy

# --- TERMINAL STYLING CONSTANTS ---
module Style
const HEADER = "\033[95m"
const BLUE = "\033[94m"
const CYAN = "\033[96m"
const GREEN = "\033[92m"
const YELLOW = "\033[93m"
const RED = "\033[91m"
const BOLD = "\033[1m"
const DIM = "\033[2m"
const END = "\033[0m"
end

"""
    to_superscript(n) -> String

Convert a non-negative integer to its unicode superscript form (e.g. 12 -> "¹²").
"""
function to_superscript(n::Integer)
    supers = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹']
    return join(supers[parse(Int, d) + 1] for d in string(n))
end

"""
    format_poly_clean(poly, var_char="x") -> String

Render a single-variable SymPy polynomial as one readable line, e.g.
"x³ - 4x + 1", instead of SymPy's default "+ -4*x" style output.
"""
function format_poly_clean(poly, var_char::String = "x")
    var = symbols(var_char)
    p = sympy.Poly(poly, var)
    term_list = p.terms()

    isempty(term_list) && return "0"

    terms = String[]
    for (power_tuple, coeff) in term_list
        power = Int(power_tuple[1])   # Sym -> native Julia Int (exponent)
        coeff_val = N(coeff)          # concrete numeric value for sign/comparisons

        # Leading term has no '+', later terms get ' + ' / ' - ' depending on sign.
        sign = if isempty(terms)
            coeff_val < 0 ? "-" : ""
        else
            coeff_val < 0 ? " - " : " + "
        end

        abs_coeff = abs(coeff_val)

        # Omit an explicit coefficient of 1 unless the term is the constant.
        coeff_str = (abs_coeff == 1 && power != 0) ? "" : string(Int(abs_coeff))

        # x^0 -> "", x^1 -> "x", x^n -> "x" * superscript(n).
        var_str = if power == 0
            ""
        elseif power == 1
            var_char
        else
            var_char * to_superscript(power)
        end

        push!(terms, sign * coeff_str * var_str)
    end

    return join(terms)
end

"""
    get_division_polynomials(n_val, x, y, A, B) -> (psi_n, phi_n)

Return the n-th division polynomial `psi_n` and the associated x-coordinate
numerator polynomial `phi_n`, for the curve y^2 = x^3 + Ax + B.

psi_n satisfies the standard recurrence (Silverman III.4):
    psi_0 = 0, psi_1 = 1, psi_2 = 2y,
    psi_3 = 3x^4 + 6Ax^2 + 12Bx - A^2,
    psi_4 = 4y(x^6 + 5Ax^4 + 20Bx^3 - 5A^2x^2 - 4ABx - 8B^2 - A^3),
and for m >= 2:
    psi_{2m}   = psi_m * (psi_{m+2} psi_{m-1}^2 - psi_{m-2} psi_{m+1}^2) / (2y)
    psi_{2m+1} = psi_{m+2} psi_m^3 - psi_{m-1} psi_{m+1}^3
phi_n is the numerator of the x-coordinate of the n-torsion multiple:
    phi_0 = 1, phi_1 = x, phi_n = x*psi_n^2 - psi_{n+1}*psi_{n-1}.

Results are memoized per call since the recurrence revisits the same index
many times (e.g. computing psi_n typically needs psi_{n-2..n+2}).
"""
function get_division_polynomials(n_val::Integer, x, y, A, B)
    psi_cache = Dict{Int, Any}()

    function psi(n::Integer)
        haskey(psi_cache, n) && return psi_cache[n]
        n == 0 && return Sym(0)
        n == 1 && return Sym(1)
        n == 2 && return 2y
        n == 3 && return 3x^4 + 6A * x^2 + 12B * x - A^2
        n == 4 && return 4y * (x^6 + 5A * x^4 + 20B * x^3 - 5A^2 * x^2 - 4A * B * x - 8B^2 - A^3)

        m = n ÷ 2
        res = if iseven(n)
            (psi(m) * (psi(m + 2) * psi(m - 1)^2 - psi(m - 2) * psi(m + 1)^2)) / (2y)
        else
            psi(m + 2) * psi(m)^3 - psi(m - 1) * psi(m + 1)^3
        end
        psi_cache[n] = expand(res)
        return res
    end

    function phi(n::Integer)
        n == 0 && return Sym(1)
        n == 1 && return x
        return expand(x * psi(n)^2 - psi(n + 1) * psi(n - 1))
    end

    return psi(n_val), phi(n_val)
end

function safe_parse_int(prompt::String)
    print(prompt)
    line = readline()
    return tryparse(Int, strip(line))
end

function main()
    println("\n\n")
    println(Style.HEADER * Style.BOLD * "╔════════════════════════════════════════════════════╗")
    println("║   ELLIPTIC CURVE DIVISION POLYNOMIAL CALCULATOR    ║")
    println("║           For Complex Multiplication               ║")
    println("╚════════════════════════════════════════════════════╝" * Style.END)

    # --- STEP 1: DEFINE CURVE ---
    println("\n" * Style.BOLD * "Step 1: Define the Curve" * Style.END)
    println("Format: y² = x³ + " * Style.CYAN * "A" * Style.END * "x + " * Style.CYAN * "B" * Style.END)

    A_in = safe_parse_int("   Enter A: " * Style.CYAN)
    print(Style.END)
    B_in = safe_parse_int("   Enter B: " * Style.CYAN)
    print(Style.END)
    (A_in === nothing || B_in === nothing) && return

    # A curve is singular (not elliptic) iff its discriminant vanishes.
    disc = 4A_in^3 + 27B_in^2
    if disc == 0
        println("\n" * Style.RED * "[CRITICAL ERROR] Singular Curve." * Style.END)
        return
    end

    cm_type = nothing
    s_symbol = ""

    # --- CM DETECTION & DESCRIPTION BLOCK ---
    # Only the two simplest CM cases are handled: j = 1728 (Z[i], B = 0)
    # and j = 0 (Z[omega], A = 0).
    if B_in == 0 && A_in != 0
        cm_type = "Gaussian"
        s_symbol = "i"
        println("\n   " * Style.GREEN * "✔ Detected CM by Gaussian Integers ℤ[i]" * Style.END)
        println("   " * Style.DIM * "----------------------------------------" * Style.END)
        println("   • " * Style.BOLD * "Ring:" * Style.END * "   ℤ[i] where i² = -1")
        println("   • " * Style.BOLD * "Map s:" * Style.END * "  Corresponds to [i]")
        println("   • " * Style.BOLD * "Action:" * Style.END * " (x, y) ↦ (-x, iy)")
        println("   " * Style.DIM * "----------------------------------------" * Style.END)

    elseif A_in == 0 && B_in != 0
        cm_type = "Eisenstein"
        s_symbol = "ω"
        println("\n   " * Style.GREEN * "✔ Detected CM by Eisenstein Integers ℤ[ω]" * Style.END)
        println("   " * Style.DIM * "----------------------------------------" * Style.END)
        println("   • " * Style.BOLD * "Ring:" * Style.END * "   ℤ[ω] where ω = (-1 + √-3)/2")
        println("   • " * Style.BOLD * "Map s:" * Style.END * "  Corresponds to [ω]")
        println("   • " * Style.BOLD * "Action:" * Style.END * " (x, y) ↦ (ζ₃x, y)")
        println("           (where ζ₃ is a primitive cube root of unity)")
        println("   " * Style.DIM * "----------------------------------------" * Style.END)

    else
        println("\n" * Style.YELLOW * "⚠ Warning: No simple CM detected." * Style.END)
        return
    end

    # --- STEP 2: DEFINE MAP ---
    println("\n" * Style.BOLD * "Step 2: Define the Map π = n + m$(s_symbol)" * Style.END)
    n_map = safe_parse_int("   Enter n (Integer): " * Style.CYAN)
    print(Style.END)
    m_map = safe_parse_int("   Enter m (Complex): " * Style.CYAN)
    print(Style.END)
    (n_map === nothing || m_map === nothing) && return

    # Computation
    x, y = symbols("x y")
    curve_eq = x^3 + A_in * x + B_in

    """Simplify and reduce modulo the curve equation y^2 = x^3 + Ax + B."""
    clean_poly(expr) = subs(simplify(expr), y^2 => curve_eq)

    println("\n" * Style.BLUE * " ➤ Computing..." * Style.END)

    local final_eq

    # --- EDGE CASE HANDLING ---
    if m_map == 0
        # Case 1: Pure Integer Map [n] -- kernel is just psi_n.
        final_poly, _ = get_division_polynomials(abs(n_map), x, y, A_in, B_in)
        final_eq = clean_poly(final_poly)

    elseif n_map == 0
        # Case 2: Pure Complex Map [m]s -- kernel is psi_m.
        final_poly, _ = get_division_polynomials(abs(m_map), x, y, A_in, B_in)
        final_eq = clean_poly(final_poly)

    else
        # Case 3: Mixed Map [n] + [m]s -- combine via the addition law.
        psi_n, phi_n = get_division_polynomials(abs(n_map), x, y, A_in, B_in)
        psi_m, phi_m = get_division_polynomials(abs(m_map), x, y, A_in, B_in)

        if cm_type == "Gaussian"
            term1 = clean_poly(phi_n * psi_m^2)
            term2 = clean_poly(phi_m * psi_n^2)
            final_eq = term1 + term2
        else # Eisenstein
            zeta = sympy.exp(2 * PI * IM / 3)
            term1 = clean_poly(phi_n * psi_m^2)
            term2 = clean_poly(phi_m * psi_n^2)
            final_eq = term1 - zeta * term2
        end
    end

    result = factor(final_eq)

    # --- OUTPUT ---
    println("\n" * "="^60)
    println(Style.BOLD * "KERNEL POLYNOMIAL Ψ_π(x)" * Style.END)
    println("(For π = $(n_map) + $(m_map)$(s_symbol))")
    println("-"^60)

    try
        clean_str = format_poly_clean(result)
        println("\n" * Style.GREEN * Style.BOLD * clean_str * Style.END * "\n")
    catch
        # format_poly_clean assumes a single-variable polynomial in x;
        # fall back to SymPy's own pretty-printer if that assumption fails.
        println(sympy.pretty(result))
    end

    println("-"^60)
    println(Style.YELLOW * "LATEX CODE:" * Style.END)
    println("\\psi_{$(n_map)+$(m_map)$(s_symbol)}(x) = " * string(sympy.latex(result)))
    println("="^60 * "\n")
end

main()
