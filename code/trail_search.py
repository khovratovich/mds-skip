"""
trail_search.py — Width-k=1 skipping trail search for MDS submatrices.

For a given t×t MDS matrix M (branch number b = t+1), this script:
  1. Extracts the (t-1)×(t-1) submatrix by deleting one row/column.
  2. Reports the theoretical minimum width from Lemma 1:  k >= (b-2)/2.
  3. Searches for arity-2 eigenvectors compatible with a k=1 skipping trail:
       - For each pair (l, r) with l < r, and each ρ with ρ^{d-1} = 1,
         check whether a = e_l + ρ*e_r satisfies M*a || a (i.e. M*a = λ*a).
  4. Reports all solutions found (or confirms absence).

Usage:
    python trail_search.py          # uses default: KoalaBear p, t=16 Cauchy, delete row/col 15
    python trail_search.py t=4      # use 4×4 Cauchy, extract 3×3 submatrix

Reference: Lemma 1 (lem:branch) and Proposition 2 (prop:mds-max2) in main.tex.
"""

import sys
import os
import itertools

# ---------------------------------------------------------------------------
# Path setup: import poseidon-tools from sibling directory
# ---------------------------------------------------------------------------
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
POSEIDON_TOOLS = os.path.join(REPO_ROOT, "..", "poseidon-tools")
sys.path.insert(0, POSEIDON_TOOLS)

from poseidon.mds_matrix import generate_mds_matrix, _mat_vec_mul

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
P = 2130706433        # KoalaBear: 2^31 - 2^24 + 1
D = 3                 # S-box exponent
T_FULL = 16           # full Poseidon1 state width
DELETE_IDX = T_FULL - 1   # row/column to delete (0-indexed)

# Override t from command line: "t=N"
for arg in sys.argv[1:]:
    if arg.startswith("t="):
        T_FULL = int(arg[2:])

# ---------------------------------------------------------------------------
# Step 1: Generate the full MDS matrix and extract the submatrix
# ---------------------------------------------------------------------------
print(f"=== Poseidon1 Cauchy MDS, t={T_FULL}, p={P}, d={D} ===\n")

M_full = generate_mds_matrix(T_FULL, P)
t = T_FULL - 1
del_i = T_FULL - 1   # delete last row and last column

M = [
    [M_full[i][j] for j in range(T_FULL) if j != del_i]
    for i in range(T_FULL) if i != del_i
]

print(f"Extracted {t}×{t} submatrix by deleting row/col {del_i}.")
print(f"Branch number of full matrix: b = {T_FULL + 1}")
print(f"Branch number of submatrix  : b = {t + 1}  (MDS, all square submatrices invertible)\n")

# ---------------------------------------------------------------------------
# Step 2: Theoretical lower bound on k from Lemma 1
# ---------------------------------------------------------------------------
b = t + 1             # MDS branch number for the t×t submatrix
k_min_theory = (b - 2 + 1) // 2   # ceiling((b-2)/2)

print(f"=== Theoretical bound (Lemma 1) ===")
print(f"  k_i >= (b-2)/2 = ({b}-2)/2 = {(b-2)/2:.1f}")
print(f"  => minimum trail width: k >= {k_min_theory}")
print(f"  => k=1 is {'IMPOSSIBLE (theory)' if k_min_theory > 1 else 'potentially possible'}\n")

# ---------------------------------------------------------------------------
# Step 3: Solve ρ^{d-1} ≡ 1 (mod p) — find valid ratios
# ---------------------------------------------------------------------------
# For d=3: ρ^2 = 1 mod p  =>  ρ ∈ {1, p-1}
# General: find all (d-1)-th roots of unity mod p.
def roots_of_unity(exp, p):
    """Return all x in GF(p) with x^exp = 1."""
    # |GF(p)^*| = p-1; the group of (exp)-th roots has order gcd(exp, p-1).
    g = pow(exp, 1, p)  # just exp mod p
    order = 1
    x = 1
    roots = []
    # Brute-force for small exp (exp = d-1 <= a few)
    for candidate in range(1, p):
        if pow(candidate, exp, p) == 1:
            roots.append(candidate)
        # Only need gcd(exp, p-1) roots; stop early
        if len(roots) == (lambda a, b: a if b == 0 else (lambda a, b: b if a % b == 0 else None)(b, a % b) or __import__('math').gcd(a,b))(exp, p-1):
            break
    return roots

import math
n_rho = math.gcd(D - 1, P - 1)   # number of (d-1)-th roots of unity
print(f"=== Valid ratios ρ with ρ^{{d-1}} = ρ^{D-1} = 1 (mod p) ===")
print(f"  gcd(d-1, p-1) = gcd({D-1}, {P-1}) = {n_rho}")

# For d=3, d-1=2: roots are ±1
rho_values = []
for rho in range(1, P):
    if pow(rho, D - 1, P) == 1:
        rho_values.append(rho)
    if len(rho_values) == n_rho:
        break

print(f"  ρ values: {[r if r <= (P-1)//2 else r-P for r in rho_values]}  (shown as signed integers)\n")

# ---------------------------------------------------------------------------
# Step 4: Search for arity-2 eigenvectors
# For each (l, r, ρ): build a = e_l + ρ*e_r, check M*a = λ*a.
# Condition: M*a has nonzeros ONLY at positions l and r.
# ---------------------------------------------------------------------------
print(f"=== Searching for arity-2 eigenvectors (k=1 trail) over {t}×{t} submatrix ===")
print(f"  Checking C({t},2) × {n_rho} = {t*(t-1)//2 * n_rho} candidates...\n")

solutions = []
for l, r in itertools.combinations(range(t), 2):
    for rho in rho_values:
        # Build test vector a = e_l + rho * e_r
        a = [0] * t
        a[l] = 1
        a[r] = rho

        # Compute Ma
        Ma = _mat_vec_mul(M, a, P)

        # Check: Ma[i] == 0 for all i not in {l, r}
        inactive_zero = all(Ma[i] == 0 for i in range(t) if i != l and i != r)
        if not inactive_zero:
            continue

        # Check: Ma is a scalar multiple of a, i.e. Ma[l]*rho == Ma[r]
        lam = Ma[l]   # proposed eigenvalue
        if (lam * rho) % P != Ma[r]:
            continue

        solutions.append((l, r, rho if rho <= (P-1)//2 else rho - P, lam))

if solutions:
    print(f"  FOUND {len(solutions)} solution(s):")
    for l, r, rho_signed, lam in solutions:
        print(f"    a = e_{l} + ({rho_signed})*e_{r},  eigenvalue λ = {lam}")
else:
    print(f"  No arity-2 eigenvectors found for k=1.")
    print(f"  This is consistent with the theoretical bound k >= {k_min_theory}.\n")

# ---------------------------------------------------------------------------
# Step 5: Minimum k for the submatrix (from Prop 1: b <= 2t/m + 2 => k >= (b-2)/2)
# For m=2 (maximum for MDS t>=3), k >= (b-2)/2.  Show exact bound.
# ---------------------------------------------------------------------------
print(f"=== Summary ===")
print(f"  Matrix:          {t}×{t} submatrix of {T_FULL}×{T_FULL} Cauchy MDS over F_p")
print(f"  Branch number:   b = {b}")
print(f"  d:               {D}")
print(f"  Lemma 1 bound:   k_i >= (b-2)/2 = {(b-2)/2:.1f}  =>  k >= {k_min_theory}")
print(f"  k=1 possible?    {'YES' if k_min_theory <= 1 else 'NO (theory rules it out)'}")
print(f"  Min trail width: k = {k_min_theory}")
print(f"  For m=2 (max by Prop 2), total active S-boxes per round: n_B = 2k >= {2*k_min_theory}")
