"""
find_trail_center.py — Find the center of the 2-round skipping trail for M=circ(2,1,1).

Trail parameters (§5 of the paper):
    Matrix   : M = circ(2,1,1) = [[2,1,1],[1,2,1],[1,1,2]]
    Dimension: t = 3
    S-box exp: d = 3
    Trail vec: a = (0, 1, -1)^T   (eigenvector of M with eigenvalue 1)
    Arity-2  : active pair (l,r) = (1,2), ratio rho = -1

Round structure (as in §2 of the paper): for each round i,
    W^(B) = M * W^(A)
    W^(C) = S(W^(B))         (coordinate-wise x^d)
    W^(A)_{next} = W^(C) + c_i

The state ratio conditions (Theorem 1, hypothesis (ii)) for a 2-round trail are:
    SR-1: W_1^(B)[2] = -W_1^(B)[1]
    SR-2: W_2^(B)[2] = -W_2^(B)[1]

Let W_1^(A) = (w0, w1, w2).  Define u0 = W_1^(B)[0] = 2w0+w1+w2.

Expanding SR-1:
    2*w0 + 3*(w1 + w2) = 0                                   ... (SR-1)

Expanding SR-2 (using SR-1 and round constant c1):
    2*u0^3 + 2*c1[0] + 3*(c1[1] + c1[2]) = 0                ... (SR-2)

Solution (Proposition 5.1 in the paper):
    u0  = cbrt( -(2*c1[0] + 3*(c1[1]+c1[2])) * inv(2) )     (unique since gcd(3,p-1)=1)
    w0  = 3*u0 * inv(4)
    w1  = free parameter
    w2  = -u0 * inv(2) - w1

Usage:
    python find_trail_center.py                              # p=11, zero constants, w1=1
    python find_trail_center.py --prime 2130706433           # KoalaBear, zero constants
    python find_trail_center.py --prime 11 --c1 1,2,3       # custom round constant
    python find_trail_center.py --prime 11 --c1 1,2,3 --c2 4,5,6 --w1 2
"""

import argparse
import sys


# ---------------------------------------------------------------------------
# Field arithmetic helpers
# ---------------------------------------------------------------------------

def inv(x: int, p: int) -> int:
    return pow(x, -1, p)


def cbrt(x: int, p: int) -> int:
    """Unique cube root mod p.  Requires gcd(3, p-1) == 1."""
    from math import gcd
    assert gcd(3, p - 1) == 1, f"gcd(3, p-1) = {gcd(3, p-1)} != 1; cube root not unique mod {p}"
    exp = pow(3, -1, p - 1)   # inverse of 3 mod (p-1)
    return pow(x, exp, p)


# ---------------------------------------------------------------------------
# Round operations (paper convention: L → S → RC)
# ---------------------------------------------------------------------------

M = [[2, 1, 1],
     [1, 2, 1],
     [1, 1, 2]]

D = 3


def mat_vec(M, v, p):
    t = len(M)
    return [sum(M[i][j] * v[j] for j in range(t)) % p for i in range(t)]


def sbox(v, d, p):
    return [pow(x, d, p) for x in v]


def add_rc(v, rc, p):
    return [(v[i] + rc[i]) % p for i in range(len(v))]


def full_round(state, rc, p):
    """One full round: Linear → S-box → Round Constant (paper convention)."""
    state = mat_vec(M, state, p)
    state = sbox(state, D, p)
    state = add_rc(state, rc, p)
    return state


# ---------------------------------------------------------------------------
# Core: find the trail center
# ---------------------------------------------------------------------------

def find_center(c1: list, p: int, w1_free: int = 1) -> dict:
    """
    Find the center W_1^(A) for the 2-round skipping trail.

    Args:
        c1       : Round constant for round 1 (list of 3 field elements mod p).
        p        : Prime field modulus.
        w1_free  : Value of the free parameter w1 (default 1).

    Returns:
        dict with keys:
            center      : W_1^(A) = [w0, w1, w2]
            u0          : W_1^(B)[0]
            W1B         : W_1^(B) = M * W_1^(A)
            W2A         : W_2^(A) = S(W_1^(B)) + c1
            W2B         : W_2^(B) = M * W_2^(A)
            sr1_ok      : bool — SR-1 satisfied
            sr2_ok      : bool — SR-2 satisfied
    """
    # Solve SR-2 for u0
    rhs = (-(2 * c1[0] + 3 * (c1[1] + c1[2]))) % p
    rhs = (rhs * inv(2, p)) % p
    u0 = cbrt(rhs, p)

    # Solve for w0 from u0 = (4/3)*w0  =>  w0 = (3/4)*u0
    w0 = (3 * u0 % p * inv(4, p)) % p

    # w2 from SR-1: w1 + w2 = -u0/2  (since 2w0+3(w1+w2)=0 and u0=(4/3)w0 => w1+w2=-u0/2)
    w1 = w1_free % p
    w2 = (-u0 * inv(2, p) - w1) % p

    center = [w0, w1, w2]

    # Propagate to verify
    W1B = mat_vec(M, center, p)
    sr1_ok = (W1B[2] == (p - W1B[1]) % p)

    W1C = sbox(W1B, D, p)
    W2A = add_rc(W1C, c1, p)
    W2B = mat_vec(M, W2A, p)
    sr2_ok = (W2B[2] == (p - W2B[1]) % p)

    return {
        "center": center,
        "u0":     u0,
        "W1B":    W1B,
        "W2A":    W2A,
        "W2B":    W2B,
        "sr1_ok": sr1_ok,
        "sr2_ok": sr2_ok,
    }


# ---------------------------------------------------------------------------
# Verification: trace the full coset flow for k inputs and confirm output
# ---------------------------------------------------------------------------

def verify_trail(center: list, c1: list, c2: list, p: int, n_samples: int = None):
    """
    Verify that the coset {center + x*a : x in F_p} is preserved (up to
    reparameterization) through 2 rounds.

    Checks that after 2 rounds every element of the coset maps to an element
    of the form (w_out[0], y, -y) for some y -- i.e. a coset of a=(0,1,-1)^T.

    Args:
        center    : W_1^(A) = [w0, w1, w2].
        c1, c2    : Round constants for rounds 1 and 2.
        p         : Prime.
        n_samples : Number of coset elements to check (default: min(p, 100)).
    """
    a = [0, 1, p - 1]   # (0, 1, -1) mod p

    if n_samples is None:
        n_samples = min(p, 100)

    failures = 0
    for x in range(n_samples):
        inp = [(center[i] + x * a[i]) % p for i in range(3)]
        out = full_round(inp, c1, p)
        out = full_round(out, c2, p)
        # Check out[2] == -out[1] mod p
        if out[2] != (p - out[1]) % p:
            failures += 1

    return failures == 0, failures


# ---------------------------------------------------------------------------
# Coset display
# ---------------------------------------------------------------------------

def print_cosets(center: list, c1: list, c2: list, p: int):
    """
    Print the full input coset and the full output coset as separate tables.

    Input coset : { center + x*a  :  x in F_p }   where a = (0,1,-1)
    Output coset: the image of the input coset after 2 full rounds.

    The output coset should be an affine coset of a=(0,1,-1)^T (i.e. every
    output element satisfies out[2] = -out[1] mod p).
    """
    a = [0, 1, p - 1]   # (0, 1, -1) mod p

    # Compute both cosets
    input_elements = []
    output_elements = []
    for x in range(p):
        inp = [(center[i] + x * a[i]) % p for i in range(3)]
        out = full_round(inp, c1, p)
        out = full_round(out, c2, p)
        input_elements.append(inp)
        output_elements.append(out)

    header = f"  {'x':>4}  {'element':^28}  ok?"
    divider = f"  {'-'*4}  {'-'*28}  ---"

    # --- Input coset ---
    print(f"\n{'='*50}")
    print(f"  Input coset  (p={p},  center = {[signed(x,p) for x in center]})")
    print(f"{'='*50}")
    print(header)
    print(divider)
    for x, inp in enumerate(input_elements):
        ok = (inp[2] == (p - inp[1]) % p)   # always True by construction; shown for symmetry
        s = f"({signed(inp[0],p):4}, {signed(inp[1],p):4}, {signed(inp[2],p):4})"
        print(f"  {x:>4}  {s:^28}  {'✓' if ok else '✗'}")

    # --- Output coset ---
    out_center = output_elements[0]
    print(f"\n{'='*50}")
    print(f"  Output coset  (c1={[signed(x,p) for x in c1]},  c2={[signed(x,p) for x in c2]})")
    print(f"  Output center = {[signed(x,p) for x in out_center]}")
    print(f"{'='*50}")
    print(header)
    print(divider)
    all_ok = True
    for x, out in enumerate(output_elements):
        ok = (out[2] == (p - out[1]) % p)
        if not ok:
            all_ok = False
        s = f"({signed(out[0],p):4}, {signed(out[1],p):4}, {signed(out[2],p):4})"
        print(f"  {x:>4}  {s:^28}  {'✓' if ok else '✗'}")
    print(f"\n  All outputs lie on coset of a=(0,1,-1): {'YES' if all_ok else 'NO'}")
    print(f"{'='*50}\n")


# ---------------------------------------------------------------------------
# Pretty printer
# ---------------------------------------------------------------------------

def signed(x, p):
    """Show x as signed integer in (-p/2, p/2]."""
    return x if x <= p // 2 else x - p


def print_result(res: dict, c1: list, c2: list, p: int, w1_free: int):
    print(f"\n{'='*60}")
    print(f"  2-Round Skipping Trail Center Finder")
    print(f"{'='*60}")
    print(f"  Prime p  : {p}")
    print(f"  d        : {D}")
    print(f"  c1       : {[signed(x, p) for x in c1]}")
    print(f"  c2       : {[signed(x, p) for x in c2]}")
    print(f"  w1 (free): {w1_free}")
    print()

    w = [signed(x, p) for x in res['center']]
    print(f"  Center W_1^(A) = ({w[0]}, {w[1]}, {w[2]})")
    print(f"  u0 = W_1^(B)[0] = {signed(res['u0'], p)}")
    print()

    b = [signed(x, p) for x in res['W1B']]
    print(f"  W_1^(B) = ({b[0]}, {b[1]}, {b[2]})")
    print(f"    SR-1: W_1^(B)[2] = -W_1^(B)[1]?  {'PASS' if res['sr1_ok'] else 'FAIL'}")

    b2 = [signed(x, p) for x in res['W2B']]
    print(f"  W_2^(B) = ({b2[0]}, {b2[1]}, {b2[2]})")
    print(f"    SR-2: W_2^(B)[2] = -W_2^(B)[1]?  {'PASS' if res['sr2_ok'] else 'FAIL'}")
    print()

    ok, failures = verify_trail(res['center'], c1, c2, p)
    n = min(p, 100)
    print(f"  Coset flow check ({n} samples): {'PASS' if ok else f'FAIL ({failures} failures)'}")
    print(f"{'='*60}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_triple(s: str, p: int) -> list:
    parts = [int(x.strip()) % p for x in s.split(',')]
    if len(parts) != 3:
        raise ValueError(f"Expected 3 comma-separated integers, got: {s!r}")
    return parts


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--prime', type=int, default=11,
                        help='Prime field modulus (default: 11)')
    parser.add_argument('--c1', type=str, default='0,0,0',
                        help='Round constant for round 1 as "a,b,c" (default: 0,0,0)')
    parser.add_argument('--c2', type=str, default='0,0,0',
                        help='Round constant for round 2 as "a,b,c" (default: 0,0,0)')
    parser.add_argument('--w1', type=int, default=1,
                        help='Free parameter w1 (default: 1)')
    args = parser.parse_args()

    p = args.prime
    c1 = parse_triple(args.c1, p)
    c2 = parse_triple(args.c2, p)

    res = find_center(c1, p, w1_free=args.w1)
    print_result(res, c1, c2, p, args.w1)

    if p <= 10000:
        print_cosets(res['center'], c1, c2, p)

    if not (res['sr1_ok'] and res['sr2_ok']):
        print("ERROR: state ratio conditions not satisfied.")
        sys.exit(1)


if __name__ == '__main__':
    main()
