"""Checks for "Borel Three-Colorings of Unit-Distance Graphs on Two Concentric Circles".

Run: python verify.py   (needs sympy and mpmath; writes verify_output.json next to this file)

1. Lemma 3.1: the palettes 3-color T_n for 4 <= n <= 500, n != 5; T_5 has no proper 3-coloring.
2. Compatibility and irrationality for every doubled pair q = a/n < 1/3, n <= 60:
   the cosine formula, the compatibility threshold, the unit criterion, and an exact
   rationality test for g (Kronecker: an algebraic integer whose conjugates are real
   and lie in [-2, 2] is 2cos of a rational multiple of pi).
3. Theorem 1.1 on every compatible pair with irrational g and n <= 60: the coloring is
   checked against the physical edge constraints cell by cell (interval intersections
   of shifted arcs, with g to 60 digits and the smallest boundary margin reported).
   Controls that must fail: q = 1/10 with k = 2, and with the s-phase 1/2 removed.
4. Proposition 5.1: a common coloring floor(3mx) mod 3 is searched for m <= 3000 on the
   same cell-by-cell check; it should be found exactly when 3 divides n.
5. Corollary 6.2: the exact identities for q = 1/10 in Q(sqrt 5).
"""
import itertools
import json
import math
import os
from fractions import Fraction as F

import mpmath as mp
import numpy as np
import sympy as sp

mp.mp.dps = 60
OUT = {}


# ---------- 1. palettes ----------
def palette(n):
    m_minus = {0: 3, 1: 2, 2: 4}[n % 3]
    m_plus = n - m_minus
    d = [1] * m_plus + [-1] * m_minus
    a = [0]
    for j in range(n - 1):
        a.append((a[-1] + d[j]) % 3)
    assert (a[-1] + d[n - 1]) % 3 == a[0]
    b = [(-a[j - 1] - a[j]) % 3 for j in range(n)]
    return a, b


def tn_edges(n):
    E = set()
    for j in range(n):
        E.add(frozenset({("A", j), ("A", (j + 1) % n)}))
        E.add(frozenset({("B", j), ("B", (j + 2) % n)}))
        E.add(frozenset({("A", j), ("B", j)}))
        E.add(frozenset({("A", j), ("B", (j + 1) % n)}))
    return [tuple(e) for e in E if len(e) == 2]


bad = []
for n in range(4, 501):
    if n == 5:
        continue
    a, b = palette(n)
    col = {("A", j): a[j] for j in range(n)}
    col.update({("B", j): b[j] for j in range(n)})
    if any(col[u] == col[v] for u, v in tn_edges(n)):
        bad.append(n)
E5 = tn_edges(5)
V5 = [("A", j) for j in range(5)] + [("B", j) for j in range(5)]
t5 = sum(1 for c in itertools.product(range(3), repeat=10)
         if all(dict(zip(V5, c))[u] != dict(zip(V5, c))[v] for u, v in E5))
OUT["lemma"] = {"n_checked": "4..500 except 5", "failures": bad,
                "T4_edges": len(tn_edges(4)), "T5_proper_3colorings": t5}
print("lemma: failures", bad, "| T_4 edges", len(tn_edges(4)), "| T_5 proper 3-colorings", t5)


# ---------- 2. compatibility and irrationality ----------
X, H = sp.symbols("X H")


def second_step(q):
    p = (2 * q) % 1
    if p > F(1, 2):
        p = 1 - p
    return p


def R(t):
    return 1 / (2 * mp.sin(mp.pi * t))


def cos_pi_g(q, p):
    r, s = R(mp.mpf(q.numerator) / q.denominator), R(mp.mpf(p.numerator) / p.denominator)
    return (r * r + s * s - 1) / (2 * r * s)


def in_family(n):
    if n % 2:
        return False
    m = n // 2
    f = sp.factorint(m)
    if m == 1:
        return False
    if len(f) == 1:
        return True  # m an odd prime power, or m = 2^f with f >= 1 (n = 2^e, e >= 2)
    return False


def g_is_rational(a, n):
    """Exact: minimal polynomial of h = 2cos(pi g) = x^3 - 3x + 1/x, x = 2cos(pi a/n)."""
    mx = sp.minimal_polynomial(2 * sp.cos(sp.pi * sp.Rational(a, n)), X)
    res = sp.resultant(sp.Poly(mx, X), sp.Poly(X * H - (X**4 - 3 * X**2 + 1), X))
    hval = mp.mpf(2) * cos_pi_g(F(a, n), second_step(F(a, n)))
    best = None
    for fac, _ in sp.factor_list(sp.Poly(res, H))[1]:
        v = abs(mp.mpf(str(sp.N(fac.as_expr().subs(H, sp.Float(str(hval), 60)), 50))))
        if best is None or v < best[0]:
            best = (v, fac)
    fac = sp.Poly(best[1], H)
    c = fac.all_coeffs()
    monic = abs(c[0]) == 1
    roots_ok = all(abs(complex(z).imag) < 1e-30 and -2 - 1e-30 <= complex(z).real <= 2 + 1e-30
                   for z in sp.Poly(fac).nroots(n=40))
    unit = abs(sp.Poly(mx, X).all_coeffs()[-1]) == 1
    return monic and roots_ok, unit, monic


c = sp.symbols("c")
q0_c = [z for z in sp.Poly(16 * c**4 - 12 * c**2 - 4 * c + 1, c).nroots(n=50)
        if z.is_real and 0.9 < z < 1]
q0 = mp.acos(mp.mpf(str(q0_c[0]))) / mp.pi
x0 = 2 * mp.mpf(str(q0_c[0]))  # the root of X^4-3X^2-2X+1 in (1,2)
assert abs(x0**4 - 3 * x0**2 - 2 * x0 + 1) < mp.mpf(10) ** -40
pairs, fam_mismatch, compat_mismatch, rational_g = [], [], [], []
for n in range(3, 61):
    for a in range(1, n):
        if math.gcd(a, n) != 1:
            continue
        q = F(a, n)
        if not q < F(1, 3):
            continue
        p = second_step(q)
        assert q < p <= F(1, 2)
        cg = cos_pi_g(q, p)
        x = 2 * mp.cos(mp.pi * a / n)
        assert abs(2 * cg - (x**3 - 3 * x + 1 / x)) < mp.mpf(10) ** -50
        compatible = abs(cg) <= 1
        if compatible != (mp.mpf(a) / n > q0):
            compat_mismatch.append((a, n))
        rec = {"q": f"{a}/{n}", "p": str(p), "n": n, "a": a, "compatible": bool(compatible)}
        if compatible:
            rat, unit, monic = g_is_rational(a, n)
            if (not unit) != in_family(n):
                fam_mismatch.append((a, n))
            if not unit and rat:
                fam_mismatch.append(("family member with rational g", a, n))
            rec.update(g_rational=bool(rat), x_unit=bool(unit), in_family=in_family(n),
                       g=float(mp.acos(cg) / mp.pi))
            if rat:
                rational_g.append((f"{a}/{n}", round(float(mp.acos(cg) / mp.pi), 9)))
        pairs.append(rec)
comp = [r for r in pairs if r["compatible"]]
OUT["pairs"] = {"doubled_pairs_q_lt_third_n_le_60": len(pairs), "compatible": len(comp),
                "x0": float(x0), "q0": float(q0), "compat_threshold_mismatches": compat_mismatch,
                "unit_criterion_mismatches": fam_mismatch, "rational_g": rational_g,
                "family_members": sum(1 for r in comp if r["in_family"]),
                "irrational_g": sum(1 for r in comp if not r["g_rational"])}
print("pairs", len(pairs), "compatible", len(comp), "x0", float(x0), "q0", float(q0),
      "| threshold mismatches", compat_mismatch, "| unit mismatches", fam_mismatch)
print("rational g:", rational_g, "| family members", OUT["pairs"]["family_members"],
      "| irrational g", OUT["pairs"]["irrational_g"])


# ---------- 3. cell-by-cell check of the construction ----------
def frac_dist(y):
    return abs(y - mp.nint(y))


def tm(x):
    return mp.mpf(x.numerator) / x.denominator if isinstance(x, F) else mp.mpf(x)


def check_coloring(L, cr, cs, q, p, half_g, s_phase):
    """r-cells [j, j+1)/L, s-cells [i - s_phase, i + 1 - s_phase)/L, colors cr[j], cs[i].
    Every pair of cells joined by some edge must get different colors.
    Returns (violations, smallest margin of a non-integral shift). Shifts in units of 1/L."""
    cr, cs = np.asarray(cr), np.asarray(cs)
    viol, margin = 0, mp.inf
    jobs = [(cr, cr, 0, 0, L * q), (cr, cr, 0, 0, -L * q)]
    jobs += [(cs, cs, -s_phase, -s_phase, L * p), (cs, cs, -s_phase, -s_phase, -L * p)]
    jobs += [(cr, cs, 0, -s_phase, L * half_g), (cr, cs, 0, -s_phase, -L * half_g)]
    for src, dst, so, do, sh in jobs:
        # [j+so, j+1+so) + sh meets [i+do, i+1+do) iff i - j lies strictly in (cst-1, cst+1)
        cst = tm(so) + tm(sh) - tm(do)
        k0 = int(mp.floor(cst - 1)) + 1
        offs = [k for k in (k0, k0 + 1) if cst - 1 < k < cst + 1]
        mg = frac_dist(cst)
        if mg > 0:
            margin = min(margin, mg)
        for k in offs:
            viol += int(np.count_nonzero(src == np.roll(dst, -(k % L))))
    return viol, margin


def build(n, a, k, s_phase=F(1, 2)):
    L = n * k
    pa, pb = palette(n)
    cr = [pa[j % n] for j in range(L)]
    cs = [pb[j % n] for j in range(L)]
    return L, cr, cs


results, worst = [], None
for r in comp:
    if r["g_rational"]:
        continue
    n, a = r["n"], r["a"]
    q, p = F(a, n), F(r["p"])
    half_g = mp.acos(cos_pi_g(q, p)) / mp.pi / 2
    inv = pow(a, -1, n)
    k = inv
    while frac_dist(k * half_g) >= mp.mpf(1) / (2 * n):
        k += n
    L, cr, cs = build(n, a, k)
    v, mg = check_coloring(L, cr, cs, q, p, half_g, F(1, 2))
    results.append({"q": r["q"], "k": k, "L": L, "violations": v, "min_margin": float(mg)})
    if worst is None or mg < worst:
        worst = mg
bad3 = [x for x in results if x["violations"]]
OUT["theorem"] = {"pairs_checked": len(results), "failures": bad3,
                  "max_L": max(x["L"] for x in results), "smallest_margin": float(worst),
                  "q_1_10": [x for x in results if x["q"] == "1/10"]}
print("theorem: checked", len(results), "failures", len(bad3), "max L", OUT["theorem"]["max_L"],
      "smallest margin", float(worst), "| q=1/10:", OUT["theorem"]["q_1_10"])

q, p = F(1, 10), F(1, 5)
half_g = mp.acos(cos_pi_g(q, p)) / mp.pi / 2
L, cr, cs = build(10, 1, 2)
c1 = check_coloring(L, cr, cs, q, p, half_g, F(1, 2))[0]
L, cr, cs = build(10, 1, 11)
c2 = check_coloring(L, cr, cs, q, p, half_g, F(0))[0]
OUT["controls"] = {"q1_10_k2_violations": c1, "q1_10_k11_no_half_phase_violations": c2}
print("controls (must be > 0):", OUT["controls"])
assert c1 > 0 and c2 > 0


# ---------- 4. common cyclic formula ----------
def common_ok(m, q, p, half_g):
    L = 3 * m
    col = [j % 3 for j in range(L)]
    return check_coloring(L, col, col, q, p, half_g, F(0))[0] == 0


prop = []
for r in comp:
    if r["g_rational"]:
        continue
    n, a = r["n"], r["a"]
    q, p = F(a, n), F(r["p"])
    half_g = mp.acos(cos_pi_g(q, p)) / mp.pi / 2
    found = None
    # rim conditions first (cheap exact test on {mq}), then the full cell check
    for m in range(1, 301):  # unfiltered: the full cell check for every small m
        if common_ok(m, q, p, half_g):
            found = m
            break
    for m in range(301, 3001) if found is None else ():
        y = (m * q) % 1
        if F(1, 3) <= y <= F(2, 3) and F(1, 3) <= (2 * y) % 1 <= F(2, 3) or \
           F(1, 3) <= y <= F(2, 3) and F(1, 3) <= (-2 * y) % 1 <= F(2, 3):
            if common_ok(m, q, p, half_g):
                found = m
                break
    prop.append((r["q"], n % 3 == 0, found))
mism = [x for x in prop if x[1] != (x[2] is not None)]
OUT["proposition"] = {"pairs": len(prop), "three_divides_n": sum(1 for x in prop if x[1]),
                      "mismatches": mism}
print("proposition: pairs", len(prop), "3|n", OUT["proposition"]["three_divides_n"], "mismatches", mism)


# ---------- 5. exact identities for q = 1/10 ----------
w = sp.sqrt(5)
t = (5 + w) / 10
r = (1 + w) / 2
s_ = sp.sqrt(t)
ok = {}
ok["R(1/10)=r"] = sp.simplify(1 / (2 * sp.sin(sp.pi / 10)) - r) == 0
ok["R(1/5)=s"] = sp.simplify((1 / (2 * sp.sin(sp.pi / 5)))**2 - t) == 0
ok["r^2+t-1=2rt"] = sp.expand(r**2 + t - 1 - 2 * r * t) == 0
T11 = sp.chebyshevt(11, X)
P = sp.Poly(sp.expand(T11 / X), X)  # T11(s) = s * P(s), P even
Pt = sp.expand(sum(cf * t**(deg // 2) for (deg,), cf in P.terms()))
ok["cos11g/s=19/25+22w/125"] = sp.simplify(Pt - (sp.Rational(19, 25) + 22 * w / 125)) == 0
ok["cos^2(11g)=1/2+6469w/31250"] = sp.simplify(sp.expand(t * Pt**2) - (sp.Rational(1, 2) + 6469 * w / 31250)) == 0
ok["cos^2(pi/10)=(5+w)/8"] = sp.simplify(sp.cos(sp.pi / 10)**2 - (5 + w) / 8) == 0
ok["difference"] = sp.simplify(sp.Rational(1, 2) + 6469 * w / 31250 - (5 + w) / 8 - (-sp.Rational(1, 8) + 10251 * w / 125000)) == 0
ok["4877/125000"] = -sp.Rational(1, 8) + 10251 * 2 / sp.Integer(125000) == sp.Rational(4877, 125000)
gnum = mp.acos(mp.sqrt((5 + mp.sqrt(5)) / 10)) / mp.pi
ok["19/110<g<21/110"] = mp.mpf(19) / 110 < gnum < mp.mpf(21) / 110
ok["minpoly 2cos(2pi g) = 5X^2-4"] = sp.minimal_polynomial(4 * t - 2, X) == 5 * X**2 - 4
OUT["corollary"] = {k: bool(v) for k, v in ok.items()}
OUT["corollary"]["g"] = float(gnum)
print("corollary:", OUT["corollary"])
assert all(ok.values())

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "verify_output.json"), "w") as fh:
    json.dump(OUT, fh, indent=1, default=str)
print("written verify_output.json")
