# Session History � Bounds on Round Skipping Attacks (SPN Hash Functions)

Each entry records the user's command verbatim and a summary of the agent's response, for reproducibility.

---

### [1] 2026-05-10 � "I want to write a cryptography paper that establishes bounds on round skipping attacks such as in 2026-306.pdf paper. That paper left an open question how far can we go with that approach. The core of my result is the theorem that is sketched in orignote.md (with a proof sketch) and formulated but not fully proven in main.tex. Make an agent that will write this paper with my assistance step by step. It should start with introduction with problem statement and related work, then necessary definitions and overview of schemes that are in scope (SPN primitives, particularly Poseidon1 and Poseidon2), then differential trails that we need, then theoretical bounds (that theorem with a proof). Then practice - we should find such a trail for concrete MDS matrices used in small field Poseidon1 (with t=16 and t=24). Finally, we will write code that does round skipping for concrete matrices and constants and trails."
**Action**: Explored workspace (main.tex, orignote.md, slides-ethproofs.tex, Rust code). Fetched 2026/306 abstract. Asked 5 clarifying questions (open question, trail width, field, sections order, venue). Created a phased plan covering Introduction -> Preliminaries -> Trails -> Bounds -> Concrete Trails -> Code.

---

### [2] 2026-05-10 � "Add to the plan the logging: every my command and the summary of your response should be recorded in the history.md for reproducibility"
**Action**: Added standing logging procedure to session plan. Established history.md format.

---

### [3] 2026-05-10 � "1. Poseidon1 has partial rounds as well, but we consider a generic structure where the number of partial rounds is unknown. After all, the number of full rounds should be bigger than the number of skipped rounds. We ignore the case when all the full rounds are skipped and we skip some of partial rounds as well."
**Action**: Recorded decision: analysis covers skipping of full rounds only (m <= R_F). Skipping into partial rounds is out of scope. This precondition to be stated in section 2 and the theorem.

---

### [4] 2026-05-10 � "2. Use matrices and constants from `poseidon-tools` folder"
**Action**: Located poseidon-tools at C:\Users\khovr\Documents\GitHub\Research\poseidon-tools. Extracted key parameters: p = 2130706433 (KoalaBear), alpha = d = 3, t = 16. MDS via Cauchy (generate_mds_matrix) or circulant (generate_circulant_mds_matrix). Round constants via Grain LFSR.

---

### [5] 2026-05-10 � "3. 3-round trails for concrete matrices will be a contribution of this paper (if we find them)."
**Action**: Added 3-round trail search as a paper contribution. Odd-m bound: b <= (4t-2)/(m-1); for m=3 gives b <= 2t-1, satisfied by MDS. So 3-round skipping is not ruled out by the bound alone.

---

### [6] 2026-05-10 � "we make the attack complexity a function of r_p. For the second question, use the circulant matrix first, and if we succeed, do also for a Cauchy matrix"
**Action**: Decision: section 5 complexity formula will be explicit in R_P. Concrete trail search order: circulant MDS first; if successful, repeat for Cauchy MDS.

---

### [7] 2026-05-10 � "use RP for the full Poseidon1 case, which is in the partial collision verifier"
**Action**: Read partial_collision_verifier.py. Confirmed: R_F = 8, R_P = 20, t = 16, p = KoalaBear, d = 3. All bounty instances use t=16 only; t=24 deferred.

---

### [8] 2026-05-10 � "yes"
**Action**: Wrote section 1 Introduction and section 1.1 Related Work in main.tex. Added \CICO and \corollary to macro.tex. Renamed "Main Results" -> two labeled sections (sec:trails, sec:bounds). Added bibliography (poseidon, poseidon2, skip2026, rebound, poseidoninitiative). Created history.md. Compiled cleanly.

---

### [9] 2026-05-10 � "Modify the truncated differential trail definition so that the compatibility condition is that B_i is compatible with A_{i+1}, not with A_i. And this compatibility means that there exists two sets of states, one being mapped of the other by the S layer, satisfying B_i and A_{i+1} respectively. This needs another definition what it means for a state set to satisfy the differential in the matrix form, it should be taken from the main theorem. Show me the result before proceeding to phase 2. DO not forget to log my questions and your responses"
**Action**: Added Definition (differential coset): set satisfying matrix differential A is {W_0 + A*X}. Added Definition (S-compatibility): diagonal B is S-compatible with A' if some coset of B maps to a coset of A' under coordinatewise x^d. Rewrote truncated differential trail: sequence {(A_i,B_i)} with B_i S-compatible with A_{i+1}. Fixed label collision eq:skipping3. Compiled cleanly. Showed result for review.

---

### [10] 2026-05-10 � "split the set satisfying definition into two"
**Action**: Split combined Definition into two separate definitions: (Def) differential coset; (Def) S-compatibility. Compiled cleanly.

---

### [11] 2026-05-10 � "Swap first two definitions. Then modify the truncated differential definition so that A maps to B by L if application of L to a satisfying set under A gives a set satisfying under B. for the history.md -- log my commands in full. Do this retroactively by taking them from the session history"
**Action**: Swapped Def 1 (differential coset) and Def 2 (truncated differential). Rewrote truncated differential condition 1: "for every differential coset W of A, its image L(W) is a differential coset of B." Rewrote history.md retroactively with full verbatim user commands for all entries. Compiled cleanly.
---

### [12] 2026-05-10 — "Add to an appropriate phase another proposition that follows from the main theorem, which says that if we skip m rounds, and next s rounds are described by function of degree D, and the CICO input condition is compatible with matrix A1 from the skipping trail, then the CICO problem is reduced to a system of equations of degree D"
**Action**: Added `\label{thm:skipping}` to the main theorem. Added Corollary (cor:cico, "CICO degree reduction") immediately after the theorem proof: if the CICO input constraint defines a differential coset of A_1, and the remaining s = R−m rounds are polynomial of degree ≤ D, then CICO reduces to a degree-D polynomial system. Proof: by Theorem thm:skipping the m-round image is a coset of B_{m+1}, parameterized by Y; the remaining rounds give degree-D equations in Y. Updated plan (Phase 4 step 8b). Compiled cleanly.

---

### [13] 2026-05-10 — "yes"
**Action**: Phase 2 — expanded §2 Preliminaries into four subsections: §2.1 SPN Permutations (formal round diagram with A/B/C state notation); §2.2 CICO Problem (definition + role in Poseidon security + bounty context); §2.3 Poseidon1 (mixed-round structure, KoalaBear bounty parameters p=2130706433/d=3/t=16/R_F=8/R_P=20, scope paragraph); §2.4 MDS Matrices and Branch Number (two definitions: branch number and MDS; submatrix characterization; Cauchy and circulant constructions). Compiled cleanly.

---

### [14] 2026-05-10 — "yes, proceed"
**Action**: Phase 3 — polished §3 Truncated Differentials and Skipping Trails. Added introductory paragraph explaining the three-layer structure of the formalism. Added Remark after Def 1 (row support vs. column space, rank). Added Remark after Def 2 characterizing the A→B condition as Im(MA)=Im(B) and noting that MDS forces equal rank. Fixed Def 5 (skipping trail): rewritten with clean grammar, labeled as \label{def:skipping-trail}, ratio condition labeled \label{eq:ratio}. Added §3.1 "A Worked Example: One-Round Skipping with Two Active Positions" (5 paragraphs: setup, ratio condition on center, propagation through S-box, skipping trail condition, summary linking to Theorem thm:skipping and Corollary cor:cico). Replaced old transition sentence with a cleaner bridge paragraph pointing to sec:bounds. Compiled cleanly.

---

### [15] 2026-05-10 — "We modify our results for MDS and other high-branch-number matrices: if matrix B has k zeros, then A cannot have fewer than t zeros, where t is a function of k and branch number b. This gives a lower bound on the number of active sboxes in a skipping trail"
**Action**: Added Lemma 1 (lem:branch, "Branch-number constraint on truncated differentials") in sec:bounds before Proposition 1: for any truncated differential (A,B) with branch number b, n_A + n_B ≥ b. Equivalently, if B has k zero diagonal entries, A has at most 2t−b−k zero rows (at least b−t+k active positions); for MDS: at least k+1 active positions. Proof: any non-zero column v of A has wt(v) ≤ n_A and wt(Mv) ≤ n_B (by Im(B)), so branch number gives n_A + n_B ≥ b. Added Corollary 2 (cor:consecutive, "Consecutive active S-boxes"): in an m-round skipping trail, n_{B_i} + n_{B_{i+1}} ≥ b for all i; total active S-boxes ≥ ⌈mb/2⌉. Proof: S-compatibility gives n_{A_{i+1}} ≤ n_{B_i}, then apply lem:branch. Rewrote proof of Proposition 1 to cite cor:consecutive cleanly (removed informal argument). Added \newtheorem{lemma}{Lemma} to macro.tex. Labeled Prop 1 as \label{prop:even-m}. Fixed the cases→plain equation for eq:skipping3. Compiled cleanly.

---

### [16] 2026-05-10 — "Generalize the definitions by replacing the diagonal B matrix with an arbitrary t×k matrix, and A is also t×k. Width k is a natural number between 1 and t. Do not change the worked example yet"
**Action**: Generalized Def 1 (differential coset): A ∈ F^{t×k}, X ∈ F^k, k called "width". Generalized Def 2 (truncated differential): both A,B ∈ F^{t×k}; condition rewritten as Im(MA) = Im(B); active positions = nonzero rows of B. Added Remark: MDS forces n_A = n_B (equal active position counts); diagonal B recovers classical form. Updated Def 3 (S-compatibility) and Def 4 (trail) to use consistent notation. Compiled cleanly.

---

### [17] 2026-05-10 — "arity-2 and skipping trail definitions"
**Action**: Added Def 5 (arity-2 matrix, \label{def:arity2}): B ∈ F^{t×k}, every column has exactly 2 nonzeros, every row at most 1 nonzero; pairs (ℓ_j, r_j), sets L={ℓ_j}, R={r_j}. Added Def 6 (skipping trail, \label{def:skipping-trail}): (m+1)-round truncated differential trail where (1) each B_i is arity-2 and (2) there exists a sequence of centers {W_i^{(A)}, W_i^{(B)}, W_i^{(C)}} propagating through the trail (W_i^{(B)}=MW_i^{(A)}, W_i^{(C)}=S(W_i^{(B)}), W_{i+1}^{(A)}=W_i^{(C)}+c_i). Ratio condition moved out of definition into Theorem 1 as explicit hypotheses. Compiled cleanly.

---

### [18] 2026-05-10 — "update worked example to t=4, circ(2,1,1,1)"
**Action**: Rewrote §3.1 as "A Worked Example: t=4, M=circ(2,1,1,1)" (\label{sec:example}). Computed circ(2,1,1,1)·(0,1,−1,0)^T = (0,1,−1,0)^T to verify A_1=B_1=A_2=B_2=(0,1,−1,0)^T is a valid 1-round trail. Verified arity-2 condition: one column, exactly 2 nonzeros at rows 1 and 2 (0-indexed). Verified ratio condition: center must satisfy W^{(B)}[1]=−W^{(B)}[0] (ratio ρ=−1) with d=3 giving (−1)^3=−1 ✓. Wrote 5 paragraphs (differential matrices, arity-2 check, ratio condition, center constraint, summary). Compiled cleanly.

---

### [19] 2026-05-10 — "Now write the full theorem and its proof"
**Action**: Wrote Theorem 1 (\label{thm:skipping}) with two explicit hypotheses: (i) matrix ratio condition (\label{eq:ratio}): ρ_{i,j}^d · A_{i+1}[ℓ_j,:] = A_{i+1}[r_j,:]; (ii) state ratio condition (\label{eq:state-ratio}): W_i^{(B)}[r_j] = ρ_{i,j} · W_i^{(B)}[ℓ_j]. Conclusion (\label{eq:skipping2}): input coset maps to {W_{m+1}^{(A)} + A_{m+1}·X'}. Full inductive proof in 4 steps per round: (1) linear layer via T_i s.t. MA_i = B_i·T_i; (2) S-box via state ratio giving v[r_j]=ρ_{i,j}·v[ℓ_j]; (3) identification via matrix ratio Im(C_i)=Im(A_{i+1}); (4) round constant shift. Compiled cleanly.

---

### [20] 2026-05-10 — "replace the lemma with a new statement using the M·A_i = B_i structure; strengthen the bound"
**Action**: Replaced Lemma 1 (lem:branch) with "Width lower bound from branch number": after reparameterization MA_i=B_i, arity-2 forces each column of B_i to have weight 2, so branch number gives wt(A_i[:,j]) ≥ b−2; A_i has at most t−b+2 zero rows. S-compatibility forces zero rows of A_{i+1} ⊇ zero rows of B_i. Conclusion: k_i ≥ (b−2)/2 for all i (\label{eq:width-bound}). Added Remark combining with ∑k_i ≤ t: b ≤ 2t/m+2; for MDS (b=t+1, t≥3): m ≤ 2t/(t-1) < 3, so m ≤ 2. Updated Corollary 2 (cor:consecutive): per-round n_{B_i} = 2k_i ≥ b−2; total ∑n_{B_i} ≥ m(b−2) (\label{eq:total-sbox}). Updated Proposition 1 (prop:even-m): b ≤ 2t/m+2 (\label{eq:skipping3}) for all m (not just even). Replaced old Proposition 2 with prop:mds-max2: MDS with t≥3 → no skipping trail of length m≥3 exists. Compiled cleanly.