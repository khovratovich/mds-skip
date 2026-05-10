// MDS test for a random N×N matrix over the prime field p = 2^31 - 2^24 + 1
// p = 2^64-2^32+1
//
// Definition used (common in block-cipher diffusion):
// A is MDS iff every k×k square submatrix (minor) is nonsingular for all k = 1..N.
//
// WARNING: This is exponential-time in N in the worst case (checks all minors).
//
// Build/run (example):
//   cargo new mds_test && cd mds_test
//   (replace src/main.rs with this file)
//   add dependency: rand = "0.8" in Cargo.toml
//   cargo run --release -- 4
//
// You can also pass an optional seed:
//   cargo run --release -- 8 12345

use rand::Rng;
use rand::{rngs::StdRng, SeedableRng};
use std::env;

const P: u64 = 18446744069414584321; // 2^64 - 2^32 + 1

#[inline]
fn mod_add(a: u64, b: u64) -> u64 {
    let s = a as u128 + b as u128;
    let s = if s >= P as u128 { s - P as u128 } else { s };
    s as u64
}

#[inline]
fn mod_sub(a: u64, b: u64) -> u64 {
    if a >= b { a - b } else { a + P - b }
}

#[inline]
fn mod_mul(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % (P as u128)) as u64
}

fn mod_pow(mut a: u64, mut e: u64) -> u64 {
    let mut r: u64 = 1;
    while e > 0 {

        if (e & 1) == 1 {
            r = mod_mul(r, a);
        }
        a = mod_mul(a, a);
        e >>= 1;
    }
    r
}

#[inline]
fn mod_inv(a: u64) -> u64 {
    // Fermat inverse since P is prime: a^(P-2) mod P
    // Caller must ensure a != 0.
    mod_pow(a, P - 2)
}

/// Returns true iff the k×k submatrix A[rows][cols] is nonsingular.
fn nonsingular_minor(matrix: &[Vec<u64>], rows: &[usize], cols: &[usize]) -> bool {
    let k = rows.len();
    debug_assert_eq!(k, cols.len());

    // Copy into a temporary k×k for elimination (simpler + usually faster than lots of indexing)
    let mut a = vec![vec![0u64; k]; k];
    for i in 0..k {
        for j in 0..k {
            a[i][j] = matrix[rows[i]][cols[j]];
        }
    }

    // Gaussian elimination to check full rank.
    // We don't need the actual determinant value, only whether it's zero.
    for col in 0..k {
        // Find pivot row at/under col with nonzero in this column.
        let mut pivot = None;
        for r in col..k {
            if a[r][col] != 0 {
                pivot = Some(r);
                break;
            }
        }
        let Some(piv) = pivot else {
            return false; // singular
        };

        // Swap into place if needed.
        if piv != col {
            a.swap(piv, col);
        }

        // Eliminate rows below.
        let inv_pivot = mod_inv(a[col][col]);
        for r in (col + 1)..k {
            if a[r][col] == 0 {
                continue;
            }
            let factor = mod_mul(a[r][col], inv_pivot);
            // row_r = row_r - factor * row_col
            a[r][col] = 0;
            for c in (col + 1)..k {
                let t = mod_mul(factor, a[col][c]);
                a[r][c] = mod_sub(a[r][c], t);
            }
        }
    }

    true
}

/// Calls `f` on every k-combination of indices 0..n-1. Early exit if f returns false.
fn for_each_combination<F>(n: usize, k: usize, mut f: F) -> bool
where
    F: FnMut(&[usize]) -> bool,
{
    let mut comb = Vec::with_capacity(k);

    fn rec<F>(start: usize, n: usize, k: usize, comb: &mut Vec<usize>, f: &mut F) -> bool
    where
        F: FnMut(&[usize]) -> bool,
    {
        if comb.len() == k {
            return f(comb);
        }
        // Need enough room to fill remaining positions
        let remaining = k - comb.len();
        for i in start..=(n - remaining) {
            comb.push(i);
            if !rec(i + 1, n, k, comb, f) {
                return false;
            }
            comb.pop();
        }
        true
    }

    rec(0, n, k, &mut comb, &mut f)
}

/// Full MDS test: checks all square minors of all sizes 1..N.
fn is_mds(matrix: &[Vec<u64>]) -> bool {
    let n = matrix.len();
    for k in 1..=n {
        let ok_rows = for_each_combination(n, k, |rows| {
            for_each_combination(n, k, |cols| nonsingular_minor(matrix, rows, cols))
        });
        if !ok_rows {
            return false;
        }
    }
    true
}

/// Generates an n×n Cauchy MDS matrix over F_P.
///
/// Uses x_i = i and y_j = n + j, so the (i,j) entry is 1/(x_i + y_j) = 1/(i + n + j).
/// All x_i are distinct, all y_j are distinct, and x_i + y_j = i + n + j >= n+1 > 0,
/// so the denominator is never zero for any n < P. Every square submatrix of a Cauchy
/// matrix is nonsingular (Cauchy determinant identity), so the result is always MDS.
fn cauchy_mds_matrix(n: usize) -> Vec<Vec<u64>> {
    let mut m = vec![vec![0u64; n]; n];
    for i in 0..n {
        for j in 0..n {
            // x_i + y_j = i + (n + j), which is in [n+1, 3n-2], never 0 mod P for practical n
            let denom = (i + n + j) as u64;
            m[i][j] = mod_inv(denom);
        }
    }
    m
}

// ── Polynomial arithmetic over F_P ──────────────────────────────────────────
// Polynomials are Vec<u64> in little-endian order: poly[i] = coeff of x^i.

fn poly_trim(p: &mut Vec<u64>) {
    while p.last() == Some(&0) {
        p.pop();
    }
}

fn poly_deg(p: &[u64]) -> Option<usize> {
    p.iter().rposition(|&c| c != 0)
}

fn poly_add_p(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut res = vec![0u64; a.len().max(b.len())];
    for (i, &c) in a.iter().enumerate() { res[i] = mod_add(res[i], c); }
    for (i, &c) in b.iter().enumerate() { res[i] = mod_add(res[i], c); }
    poly_trim(&mut res);
    res
}

fn poly_sub_p(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut res = vec![0u64; a.len().max(b.len())];
    for (i, &c) in a.iter().enumerate() { res[i] = c; }
    for (i, &c) in b.iter().enumerate() { res[i] = mod_sub(res[i], c); }
    poly_trim(&mut res);
    res
}

fn poly_mul_p(a: &[u64], b: &[u64]) -> Vec<u64> {
    if a.is_empty() || b.is_empty() { return vec![]; }
    let mut res = vec![0u64; a.len() + b.len() - 1];
    for (i, &ai) in a.iter().enumerate() {
        if ai == 0 { continue; }
        for (j, &bj) in b.iter().enumerate() {
            res[i + j] = mod_add(res[i + j], mod_mul(ai, bj));
        }
    }
    poly_trim(&mut res);
    res
}

/// Polynomial division: returns (quotient, remainder) such that a = q·b + r.
fn poly_divrem_p(a: &[u64], b: &[u64]) -> (Vec<u64>, Vec<u64>) {
    let db = poly_deg(b).expect("poly_divrem_p: divisor is zero");
    let inv_lead_b = mod_inv(b[db]);
    let mut rem = a.to_vec();
    poly_trim(&mut rem);
    let da = match poly_deg(&rem) {
        None => return (vec![], vec![]),
        Some(d) => d,
    };
    if da < db { return (vec![], rem); }
    let mut quot = vec![0u64; da - db + 1];
    for i in (0..=(da - db)).rev() {
        let top = i + db;
        let lead = rem.get(top).copied().unwrap_or(0);
        if lead == 0 { continue; }
        let q = mod_mul(lead, inv_lead_b);
        quot[i] = q;
        for j in 0..=db {
            let v = mod_mul(q, b[j]);
            rem[i + j] = mod_sub(rem[i + j], v);
        }
    }
    poly_trim(&mut rem);
    poly_trim(&mut quot);
    (quot, rem)
}

fn poly_rem_p(a: &[u64], b: &[u64]) -> Vec<u64> {
    poly_divrem_p(a, b).1
}

/// Polynomial GCD, returned monic (or empty for the zero polynomial).
fn poly_gcd_p(a: &[u64], b: &[u64]) -> Vec<u64> {
    if b.is_empty() {
        let mut a = a.to_vec();
        poly_trim(&mut a);
        if let Some(d) = poly_deg(&a) {
            let inv = mod_inv(a[d]);
            for c in &mut a { *c = mod_mul(*c, inv); }
        }
        return a;
    }
    let r = poly_rem_p(a, b);
    poly_gcd_p(b, &r)
}

/// Compute base^exp mod modulus in F_P[x] using repeated squaring.
fn poly_powmod_p(base: &[u64], mut exp: u64, modulus: &[u64]) -> Vec<u64> {
    let mut result = vec![1u64];
    let mut b = base.to_vec();
    while exp > 0 {
        if exp & 1 == 1 {
            result = poly_rem_p(&poly_mul_p(&result, &b), modulus);
        }
        b = poly_rem_p(&poly_mul_p(&b, &b), modulus);
        exp >>= 1;
    }
    result
}

/// Find one root of `f` in F_P.
///
/// Algorithm:
///   1. Compute h = gcd(f, x^P - x mod f), which isolates all distinct linear factors.
///   2. Recursively split h via Cantor–Zassenhaus: pick random t, compute
///      gcd(h, (x+t)^((P-1)/2) - 1) to obtain a non-trivial factor with probability ~1/2.
///
/// Returns `None` if f has no roots in F_P (degree 0 or no linear factors).
fn find_root_fp(f: &[u64], rng: &mut impl Rng) -> Option<u64> {
    let mut f = f.to_vec();
    poly_trim(&mut f);
    let d = poly_deg(&f)?;
    if d == 0 { return None; }

    // Make monic
    let inv_lead = mod_inv(f[d]);
    for c in &mut f { *c = mod_mul(*c, inv_lead); }

    // h = gcd(f, x^P - x): product of all distinct linear factors of f over F_P
    let x_poly = vec![0u64, 1u64];
    let xp = poly_powmod_p(&x_poly, P, &f);
    let xp_minus_x = poly_sub_p(&xp, &x_poly);
    let h = if xp_minus_x.is_empty() {
        f.clone() // f already divides x^P - x
    } else {
        poly_gcd_p(&f, &xp_minus_x)
    };

    let dh = poly_deg(&h)?;
    if dh == 0 { return None; } // no roots in F_P

    root_split(&h, rng)
}

/// Cantor–Zassenhaus splitting: finds one root of a square-free polynomial
/// all of whose roots lie in F_P.
fn root_split(f: &[u64], rng: &mut impl Rng) -> Option<u64> {
    let d = poly_deg(f)?;
    if d == 0 { return None; }

    // Base case: linear factor f = x + c (monic) → root = −c
    if d == 1 {
        return Some(mod_sub(0, f[0]));
    }

    let half_order = (P - 1) / 2;
    loop {
        let t: u64 = rng.gen_range(0..P);
        // g = gcd(f, (x+t)^((P−1)/2) − 1)
        // Each non-zero element of F_P satisfies a^((P-1)/2) = ±1, so this
        // splits the roots roughly in half on average.
        let x_plus_t = vec![t, 1u64];
        let powered = poly_powmod_p(&x_plus_t, half_order, f);
        let powered_m1 = poly_sub_p(&powered, &[1u64]);
        if powered_m1.is_empty() { continue; }

        let g = poly_gcd_p(f, &powered_m1);
        let dg = poly_deg(&g).unwrap_or(0);
        if dg == 0 || dg == d { continue; } // trivial split, retry

        // Recurse on both halves; return whichever yields a root first
        if let Some(r) = root_split(&g, rng) { return Some(r); }
        let (fq, _) = poly_divrem_p(f, &g);
        return root_split(&fq, rng);
    }
}

/// Scale every coefficient of a polynomial by a field element.
fn poly_scale(p: &[u64], s: u64) -> Vec<u64> {
    let mut r: Vec<u64> = p.iter().map(|&c| mod_mul(c, s)).collect();
    poly_trim(&mut r);
    r
}

/// Evaluate polynomial at `z` using Horner's method.
fn poly_eval(p: &[u64], z: u64) -> u64 {
    let mut r = 0u64;
    for &c in p.iter().rev() {
        r = mod_add(mod_mul(r, z), c);
    }
    r
}

/// Formal derivative of a polynomial over F_P: (sum a_i x^i)' = sum i*a_i x^(i-1).
fn poly_derivative(p: &[u64]) -> Vec<u64> {
    if p.len() <= 1 { return vec![]; }
    let mut d: Vec<u64> = (1..p.len())
        .map(|i| mod_mul(i as u64, p[i]))
        .collect();
    poly_trim(&mut d);
    d
}

/// Build Q(z) = ∏_i (z + x_i) from a slice of field elements.
fn build_product_poly(x: &[u64]) -> Vec<u64> {
    let mut q = vec![1u64];
    for &xi in x {
        q = poly_mul_p(&q, &[xi, 1u64]);
    }
    q
}

/// Build the polynomial p(z) = λ·Q(z) − Q'(z),
/// where Q(z) = ∏_i (z + x_i).
///
/// Its roots are exactly the y_j values satisfying:
///   ∑_i 1/(x_i + y_j) = λ   (column-j sum of the Cauchy matrix equals λ)
///
/// Derivation: Q'(z)/Q(z) = ∑_i 1/(z + x_i), so the condition Q'(z) = λ Q(z)
/// is equivalent to ∑_i 1/(z + x_i) = λ.  The polynomial is λQ − Q', degree n.
fn eigenvalue_poly(x: &[u64], lambda: u64) -> Vec<u64> {
    let q  = build_product_poly(x);
    let dq = poly_derivative(&q);
    poly_sub_p(&poly_scale(&q, lambda), &dq)
}

/// Find ALL roots of `poly` in F_P by repeated root-finding + deflation.
/// After finding each root r, divides the polynomial by (z − r) before searching again.
fn find_all_roots(poly: &[u64], rng: &mut impl Rng) -> Vec<u64> {
    let mut roots = Vec::new();
    let mut p = poly.to_vec();
    poly_trim(&mut p);
    loop {
        match poly_deg(&p) {
            None | Some(0) => break,
            _ => {}
        }
        match find_root_fp(&p, rng) {
            None => break, // no more roots in F_P
            Some(r) => {
                roots.push(r);
                // Deflate: p ← p / (z − r)
                let (q, _rem) = poly_divrem_p(&p, &[mod_sub(0, r), 1u64]);
                p = q;
            }
        }
    }
    roots
}

/// Generates an n×n Cauchy MDS matrix with a given vector eigenvalue.
///
/// `eigenvalues[j]` is the required sum of column `j`: ∑_i 1/(x_i + y_j) = λ_j.
///
/// Algorithm:
///   1. Draw n distinct x_i at random.
///   2. For each j, build p_j(z) = λ_j·Q(z) − Q'(z)  where Q = ∏_i(z+x_i).
///      The root of p_j is the unique y_j satisfying the column-j sum condition.
///   3. Find one valid root y_j per column (Cauchy condition x_i+y_j≠0, distinct y_j).
///   4. Verify MDS; retry from step 1 on failure.
///
/// Panics if `eigenvalues.len() != n`.
/// Returns `(matrix, x, y)` where `x` and `y` are the Cauchy parameters.
fn cauchy_mds_with_eigenvalue(n: usize, eigenvalues: &[u64], rng: &mut impl Rng) -> (Vec<Vec<u64>>, Vec<u64>, Vec<u64>) {
    assert_eq!(eigenvalues.len(), n, "eigenvalues length must equal n");
    use std::collections::HashSet;
    loop {
        // ── Step 1: generate n distinct nonzero x_i ──────────────────────────
        let mut x: Vec<u64> = Vec::with_capacity(n);
        let mut x_set: HashSet<u64> = HashSet::new();
        while x.len() < n {
            let xi: u64 = rng.gen_range(1..P);
            if x_set.insert(xi) { x.push(xi); }
        }

        // Precompute Q(z) and Q'(z) once; they are the same for every column.
        let q  = build_product_poly(&x);
        let dq = poly_derivative(&q);

        // ── Step 2 & 3: for each column j, solve p_j(y_j) = 0 ───────────────
        let mut y: Vec<u64> = Vec::with_capacity(n);
        let mut y_set: HashSet<u64> = HashSet::new();
        let mut ok = true;

        'col: for j in 0..n {
            let lj = eigenvalues[j];
            // p_j(z) = λ_j · Q(z) − Q'(z)
            let pj = poly_sub_p(&poly_scale(&q, lj), &dq);
            if poly_deg(&pj).map_or(true, |d| d == 0) { ok = false; break; }

            // Collect all F_P-roots of p_j, then pick the first valid one.
            let candidates = find_all_roots(&pj, rng);
            for yj in candidates {
                if !y_set.contains(&yj)
                    && x.iter().all(|&xi| mod_add(xi, yj) != 0)
                    // sanity: column sum must equal λ_j
                    && x.iter().fold(0u64, |acc, &xi| mod_add(acc, mod_inv(mod_add(xi, yj)))) == lj
                {
                    y_set.insert(yj);
                    y.push(yj);
                    continue 'col;
                }
            }
            // No valid y_j found for this column with these x_i — retry.
            ok = false;
            break;
        }

        if !ok || y.len() < n { continue; }

        // ── Step 4: build the matrix and verify MDS ───────────────────────────
        let mut m = vec![vec![0u64; n]; n];
        for i in 0..n {
            for j in 0..n {
                m[i][j] = mod_inv(mod_add(x[i], y[j]));
            }
        }

        if is_mds(&m) {
            return (m, x, y);
        }
    }
}

fn random_matrix(n: usize, rng: &mut impl Rng) -> Vec<Vec<u64>> {
    let mut m = vec![vec![0u64; n]; n];
    for i in 0..n {
        for j in 0..n {
            m[i][j] = rng.gen_range(0..64);
        }
    }
    m
}

fn random_circulant_matrix(n: usize, rng: &mut impl Rng) -> Vec<Vec<u64>> {
    let mut first_row = vec![0u64; n];
    for j in 0..n {
        first_row[j] = rng.gen_range(1..8); // Avoid zero to increase chance of MDS (not guaranteed)
    }
    let mut m = vec![vec![0u64; n]; n];
    for i in 0..n {
        for j in 0..n {
            m[i][j] = first_row[(j + i) % n];
        }
    }
    m
}

#[test]
fn test_matrix() -> Result<(), ()> {
    let mut first_row = vec![1, 1, 2, 1, 8, 9, 10, 7, 5, 9, 4, 10];
    let n = first_row.len();
    let mut m = vec![vec![0u64; n]; n];
    for i in 0..n {
        for j in 0..n {
            m[i][j] = first_row[(j + i) % n];
        }
    }

    let ok = is_mds(&m);
        if ok {
            println!("Test success!");
            // Print a small matrix (optional; comment out if you want speed only)
            if n <= 16 {
                println!("Matrix (mod P={}):", P);
                for row in &m {
                    println!("{:?}", row);
                }
                println!();
            }
            Ok(())
        }
        else {
            Err(())
        }
}

#[test]
fn test_find_root_fp() {
    let mut rng: StdRng = StdRng::seed_from_u64(42);

    // (x - 7)(x - 13) = x^2 - 20x + 91
    let r1: u64 = 7;
    let r1_neg = mod_sub(0, r1);
    let r2: u64 = 13;
    let r2_neg = mod_sub(0, r2);
    // f = (x + (-7))(x + (-13)) = x^2 + (-20)x + 91
    let neg_sum = mod_add(r1_neg, r2_neg);   // -(7+13) = -20
    let product = mod_mul(r1_neg, r2_neg);    // (-7)(-13) = 91
    let f = vec![product, neg_sum, 1u64];     // 91 - 20x + x^2
    let root = find_root_fp(&f, &mut rng).expect("should find a root");
    assert!(root == 7 || root == 13, "unexpected root {root}");
    println!("quadratic root: {root} ✓");

    // Degree-5 polynomial with known root 42: f = (x - 42)(x^4 + 1)
    // x^4 + 1 has no roots in F_P (check: a^4 = -1 has no solution since P ≡ 1 mod 8 ...
    // but let's just test that we find 42)
    // (x - 42)(x + 0)(x + 1)(x + 2)(x + 3) — five known roots
    let roots = [42u64, 100, 200, 300, 400];
    let mut poly = vec![1u64]; // start with 1
    for &r in &roots {
        // multiply by (x - r) = (x + (-r))
        let neg_r = mod_sub(0, r);
        poly = poly_mul_p(&poly, &[neg_r, 1u64]);
    }
    let found = find_root_fp(&poly, &mut rng).expect("should find a root");
    assert!(roots.contains(&found), "unexpected root {found}");
    println!("degree-5 root: {found} ✓");

    // No-root case: f = x^2 + 1 (irreducible if -1 is not a QR mod P)
    // P ≡ 1 mod 4 so -1 IS a QR; skip no-root test and just test a constant
    let no_root_poly = vec![5u64]; // constant 5, no roots
    assert!(find_root_fp(&no_root_poly, &mut rng).is_none());
    println!("constant polynomial: no root ✓");
}

#[test]
fn test_cauchy_mds() {
    for n in 1..=8 {
        let m = cauchy_mds_matrix(n);
        assert!(is_mds(&m), "cauchy_mds_matrix({n}) is not MDS");
        println!("cauchy_mds_matrix({n}): MDS ✓");
    }
}

#[test]
fn test_cauchy_mds_with_eigenvalue() {
    let mut rng: StdRng = StdRng::seed_from_u64(99);

    // Eigenvector (1, 1, -1, -1, 0, 0, 0) in F_P:
    //   col sums = [1, 1, P-1, P-1, 0, 0, 0]
    let neg1 = mod_sub(0, 1);
    let eigenvalues: Vec<u64> = vec![1, 1, neg1, neg1, 0, 0, 0];
    let n = eigenvalues.len(); // 7

    let (m, x, y) = cauchy_mds_with_eigenvalue(n, &eigenvalues, &mut rng);

    // Verify MDS
    assert!(is_mds(&m), "eigenvector (1,1,-1,-1,0,0,0): matrix is not MDS");

    // Verify each column sum matches the requested eigenvector entry
    for j in 0..n {
        let col_sum = (0..n).fold(0u64, |acc, i| mod_add(acc, m[i][j]));
        let expected = eigenvalues[j];
        assert_eq!(col_sum, expected, "col {j}: sum={col_sum} ≠ expected={expected}");
    }

    // Pretty-print with signed integers
    let signed = |v: u64| -> i128 {
        if v <= P / 2 { v as i128 } else { v as i128 - P as i128 }
    };
    println!("Cauchy MDS (7×7) with eigenvector (1,1,-1,-1,0,0,0): MDS ✓");
    println!("x_i: {:?}", x.iter().map(|&v| signed(v)).collect::<Vec<_>>());
    println!("y_j: {:?}", y.iter().map(|&v| signed(v)).collect::<Vec<_>>());
    println!("Column sums: {:?}", eigenvalues.iter().map(|&v| signed(v)).collect::<Vec<_>>());
    println!("Matrix (signed):");
    for row in &m {
        println!("  {:?}", row.iter().map(|&v| signed(v)).collect::<Vec<_>>());
    }
}

fn random_small_weight_circulant_matrix(n: usize, rng: &mut impl Rng) -> Vec<Vec<u64>> {
    //generate a vector of powers of two up to 2^(w-1) where w is the desired weight (number of nonzero entries in the first row)
    let w = 5; // Desired weight (number of nonzero entries in the first row)
    let mut powers_of_two = Vec::with_capacity(w);
    for i in 0..w {
        powers_of_two.push(1 << i);
    }
    //for each position in the first row, choose a random sum of k distinct powers of two, where k is a random number between 1 and w (inclusive)
    let k = 2; // Max number of nonzero entries in the first row
    let mut first_row = vec![0u64; n];
    for j in 0..n {
        for _ in 0..rng.gen_range(1..=k) {
            let pow = powers_of_two[rng.gen_range(0..w)];
            first_row[j] = mod_add(first_row[j], pow);
        }
    }
    let mut m = vec![vec![0u64; n]; n];
    for i in 0..n {
        for j in 0..n {
            m[i][j] = first_row[(j + i) % n];
        }
    }
    m
}


fn main() {


    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} N [seed]", args[0]);
        eprintln!("Example: {} 4", args[0]);
        std::process::exit(2);
    }

    let n: usize = args[1].parse().expect("N must be an integer");
    if n == 0 {
        eprintln!("N must be >= 1");
        std::process::exit(2);
    }

    let mut rng: StdRng = if args.len() >= 3 {
        let seed: u64 = args[2].parse().expect("seed must be an integer");
        StdRng::seed_from_u64(seed)
    } else {
        // Use an OS-seeded RNG to create a deterministic StdRng seed.
        let seed: u64 = rand::random();
        StdRng::seed_from_u64(seed)
    };

    for i in 0..10000 {
        let m = random_small_weight_circulant_matrix(n, &mut rng);
        let ok = is_mds(&m);
        if ok {
            println!("Found MDS matrix on iteration {}!", i);
            // Print a small matrix (optional; comment out if you want speed only)
            if n <= 16 {
                println!("Matrix (mod P={}):", P);
                for row in &m {
                    println!("{:?}", row);
                }
                println!();
            }
            break;
        }
    }



    //println!("N={}  =>  MDS: {}", n, ok);
}