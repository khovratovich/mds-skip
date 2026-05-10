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