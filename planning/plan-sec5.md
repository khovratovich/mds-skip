# Plan: §5 Concrete 2-Round Skipping Trail

**Date**: 2026-05-12
**Scope**: §5 of tex/main.tex + Introduction contributions list fix.
**Decisions**: F_11 as primary example, then KoalaBear; t=3, k=1; zero round constants; t=16 deferred.

---

## A. Fix Introduction (contributions list)

Second bullet (tight bounds): replace stale two-formula bound with b <= 2t/m+2 for
all m; update "rules out m >= 4" -> "rules out m >= 3" for MDS t >= 3; remove
"three-round skipping is not excluded" (it now IS excluded by Proposition 2).

Third bullet (concrete trails): replace the 3-round search language with:
"We exhibit an explicit tight m=2 skipping trail for circ(2,1,1) over F_11 and
over KoalaBear, showing the bound m <= 2 of Proposition 2 is achieved."

---

## B. §5 Structure

  \section{Concrete Trails}\label{sec:concrete}
    opening paragraph
    \subsection{The MDS Matrix}
    \subsection{The Skipping Trail}
    \subsection{Center Constraints and Coset Flow}
    \subsection{Extension to KoalaBear}

---

## §5.1 — The MDS Matrix

M = circ(2,1,1) = [[2,1,1],[1,2,1],[1,1,2]]

MDS verification over F_11:
  det(M) = 6 - 1 - 1 = 4 != 0 mod 11.
  All 2x2 minors: det[[2,1],[1,2]]=3, det[[1,1],[2,1]]=-1, all nonzero mod 11.
  Key: c0^2 - c1^2 = 3 != 0 mod 11.
  Branch number b = t+1 = 4.  M is MDS.

---

## §5.2 — The Trail

a = (0, 1, -1)^T.

Eigenvector: M*a = (0+1-1, 0+2-1, 0+1-2)^T = (0,1,-1)^T = a.  Eigenvalue 1.
Arity-2: 2 nonzeros at rows 1,2; each row at most 1 nonzero.
Active pair (l,r) = (1,2), ratio rho = -1.
Matrix ratio condition: rho^d * a[l] = (-1)^3 * 1 = -1 = a[r].  Holds for all odd d.

Trail: A_i = B_i = a for i = 1, 2, 3.

---

## §5.3 — Center Constraints and Coset Flow (zero round constants)

State ratio condition for round i: W_i^(B)[2] = -W_i^(B)[1].

Center: W_1^(A) = (0, 1, -1)^T over F_11 (c=1).

Coset flow:
  Input:     {(0, 1+x, -1-x) : x in F_11}
  After L_1: {(0, 1+x, -1-x)}              [M fixes eigenspace, eigenvalue 1]
  After S_1: {(0, (1+x)^3, -(1+x)^3)}     [d=3 odd preserves sign flip]
  After L_2: {(0, (1+x)^3, -(1+x)^3)}
  After S_2: {(0, (1+x)^9, -(1+x)^9)}     [bijection: gcd(9,10)=1 over F_11]

Output coset is a coset of a centered at (0,1,-1)^T.  2-round skipping trail confirmed.

Proposition: For M=circ(2,1,1), d=3, zero round constants, any c != 0:
  center (0,c,-c)^T gives a valid 2-round skipping trail of width k=1.

Remark: Proposition 2 (prop:mds-max2) shows m >= 3 is impossible. So m=2 is tight.

---

## §5.4 — Extension to KoalaBear

p = 2^31 - 2^24 + 1 = 2130706433.

  det(M) = 4 != 0 mod p.                       PASS.
  3 != 0 mod p.                                 PASS.  => M is MDS.
  gcd(3, p-1): p-1 = 2^24 * 127, 127 = 42*3+1, so gcd(3,p-1)=1.  PASS.
  => x |-> x^3 is a bijection on F_p.

Same trail a=(0,1,-1)^T works over KoalaBear (zero round constants).
Non-zero constants: center still exists (unique cube root), deferred to t=16 section.

---

## Verification

1. pdflatex -interaction=nonstopmode main.tex  -- no errors.
2. M*(0,1,-1)^T = (0+1-1, 0+2-1, 0+1-2)^T = (0,1,-1)^T.  PASS.
3. (-1)^3 = -1 in any field.  PASS.
4. gcd(9,10)=1 so x^9 bijection on F_11.  PASS.
5. det=4, 3 != 0 mod 11.  PASS.

---

## Files to Modify

- tex/main.tex: contributions list (intro) + §5 before \begin{thebibliography}.
- tex/macro.tex: no changes.

---

## Out of Scope (deferred)

- Poseidon1-KoalaBear t=16 trails.
- Non-zero round constants for t=3.
- k >= 2 trails.
- Rust implementation.