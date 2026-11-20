# Borel Three-Colorings of Unit-Distance Graphs on Two Concentric Circles

Yitzchak Shmalo

The paper and the script that checks its computations.

Take two concentric circles whose radii make unit chords subtend the rational numbers of turns q < p, and join two points of their union when they are at distance one. Suppose that p = ±2q (mod 1) and that the angle πg at the centre between the endpoints of a unit segment joining the two circles is an irrational multiple of π. Then the chromatic number, the Lebesgue measurable chromatic number and the Borel chromatic number of this graph all equal 3. The coloring is constant on finitely many arcs and comes from a 3-coloring of one of Wilson's rose window graphs. Under the same hypotheses, the coloring x ↦ ⌊3mx⌋ mod 3 of the turn x on both circles is proper for some m exactly when 3 divides the denominator of q. The angle g is irrational for every compatible pair with p = ±2q except q = 1/5, 1/7, 2/7, 1/9, 2/9, 1/12, so the result holds for all other pairs.

## Contents

- `paper/main.tex`, `paper/main.pdf`: the paper.
- `code/verify.py`: checks Lemma 3.1 for 4 ≤ n ≤ 500, the compatibility and irrationality of g for every pair q = a/n < 1/3 with n ≤ 60, the coloring of Theorem 1.1 on every such pair with irrational g, Proposition 5.1, and the identities of Corollary 6.2.
- `code/verify_output.json`: its output.
- `code/rational_g.py`: the computation in the proof of Proposition 6.3 for denominators above 60. It applies Mann's theorem to the 41 partitions of the six terms, which shows that the denominator divides 15 times one of 2, 4, 6, 10, 14, 18. For the 95 pairs left above 60 it shows that some conjugate of 2cos(πg) lies outside [-2, 2], so g is irrational. Pairs with denominator at most 60 are handled exactly in `verify.py`.
- `code/rational_g_output.json`: its output.

Requires Python 3, SymPy and mpmath. Run `python code/verify.py` and `python code/rational_g.py`.
