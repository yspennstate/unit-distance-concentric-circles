"""Checks for the proposition on the rational values of g (Section 6 and Computations, item 3).

Part 1. If g is rational, the six numbers lambda_1..lambda_6 below give a vanishing sum of roots of unity
sum_j exp(pi i lambda_j) = 0. Split the indices into classes of size at least two, each a minimal
vanishing sum. By Mann's theorem two terms in one class differ by a 30th root of unity, so
delta1*q + delta2*g lies in Z/15 for every (delta1, delta2) in the lattice Lambda spanned by the coefficient
differences inside the classes. For each of the 41 partitions this prints the least nu > 0 with
(nu, 0) in Lambda, so that n divides 15 nu.

Part 2. For every compatible q = a/n (q0 <= q < 1/3) with n > 60 dividing some 15 nu, it computes the
conjugates kappa_b = x_b^3 - 3 x_b + 1/x_b, x_b = 2 cos(pi a b / n), b coprime to 2n, of kappa = 2 cos(pi g).
If g were rational every conjugate would lie in [-2, 2]; the script prints the least excess over 2
of the largest conjugate. The pairs with n <= 60 are handled exactly in verify.py.
"""
import json
import math
from itertools import combinations

import mpmath as mp

mp.mp.dps = 50
# lambda_j = (coefficient of q, coefficient of g) + an integer constant
T = {1: (2, 1), 2: (-2, -1), 3: (2, -1), 4: (-2, 1), 5: (5, 0), 6: (-5, 0)}


def partitions(s):
    if not s:
        yield []
        return
    first, rest = s[0], s[1:]
    for k in range(0, len(rest) + 1):
        for c in combinations(rest, k):
            left = [x for x in rest if x not in c]
            for p in partitions(left):
                yield [(first,) + c] + p


def least_c(vecs):
    """Least c > 0 with (c, 0) in the Z-span of vecs (vecs span a lattice containing such a vector)."""
    det = 0
    for (a1, b1), (a2, b2) in combinations(vecs, 2):
        det = math.gcd(det, a1 * b2 - a2 * b1)
    ga = gb = 0
    for a, b in vecs:
        ga, gb = math.gcd(ga, a), math.gcd(gb, b)
    if gb == 0:          # every difference has beta = 0: the span is ga*Z x {0}
        return ga
    assert det > 0, "rank one with beta != 0 would give no bound on q"
    return det // gb     # index det = c * gb


parts = [P for P in partitions([1, 2, 3, 4, 5, 6]) if all(len(b) >= 2 for b in P)]
assert len(parts) == 41
cs = {}
for P in parts:
    vecs = [(T[b[0]][0] - T[j][0], T[b[0]][1] - T[j][1]) for b in P for j in b[1:]]
    c = least_c(vecs)
    assert c > 0
    cs.setdefault(c, []).append(P)
assert sorted(cs) == [2, 4, 6, 10, 14, 18]
assert [len(cs[c]) for c in sorted(cs)] == [34, 2, 1, 1, 1, 2]
Ns = sorted(15 * c for c in cs)
print("partitions:", len(parts), " n divides one of", Ns)
for c in sorted(cs):
    print("  c =", c, ":", len(cs[c]), "partitions")

q0 = mp.mpf("0.07790844892243048")   # Lemma 2.1; checked below: no a/n with n among the denominators above 60 is within 1e-4 of q0
dens = sorted({d for N in Ns for d in range(61, N + 1) if N % d == 0})
gap = min(abs(mp.mpf(a) / n - q0) for n in dens for a in range(1, n))
assert gap > mp.mpf(10) ** -4      # 7/90 is the nearest, at 1.3e-4
print("least distance from a/n (n in the denominators below) to q0:", mp.nstr(gap, 3))
rows = []
for n in dens:
    ks = [k for k in range(1, 2 * n) if math.gcd(k, 2 * n) == 1]
    for a in range(1, n):
        if math.gcd(a, n) != 1 or not (q0 <= mp.mpf(a) / n < mp.mpf(1) / 3):
            continue
        hmax = max(abs(x**3 - 3 * x + 1 / x) for x in (2 * mp.cos(mp.pi * a * k / n) for k in ks))
        rows.append((a, n, float(hmax - 2)))
excess = min(r[2] for r in rows)
print("denominators above 60:", dens)
print("pairs:", len(rows), " least excess of the largest conjugate over 2:", round(excess, 4))
assert excess > 3.04 and len(rows) == 95
json.dump({"bound_denominators": Ns, "c_counts": {str(c): len(v) for c, v in cs.items()},
           "denominators_above_60": dens, "pairs_above_60": len(rows), "least_excess": excess,
           "rows": rows}, open("rational_g_output.json", "w"), indent=1)
