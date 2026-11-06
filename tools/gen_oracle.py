#!/usr/bin/env python3
"""Generate tests/oracle_cases.hpp: scipy reference cases for expm, care, dare, the
factorizations (LU, QR, Cholesky, LDLᵀ, SVD), eigenvalues, logm, sqrtm, lyap / dlyap,
ZOH / Tustin discretization, the SISO conversions and pole placement.

Run by hand (py -3 tools/gen_oracle.py on Windows); the header is committed and
make never calls Python. Seeds are fixed, so a re-run on the same numpy/scipy
reproduces the header byte for byte.

Each case is a struct of constexpr std::array<double, ...> literals (row-major,
printed with repr, which round-trips a double exactly) plus the sizes, a
conditioning estimate for scaling the test bounds, and a float_ok flag. The
header does not include Damp, so cases are usable inside static_assert.

The scipy answer is the secondary check. The tests' primary checks are the
Riccati residual, symmetry, definiteness and closed-loop stability (G2.4), and for
the later families the residuals and reconstructions a backward-stable method meets.
Where a reference must be exact (solutions, determinants, unpivoted LDLᵀ, Leverrier
coefficients, polynomial expansions, Markov parameters) it is computed in rational
arithmetic on the stored doubles and rounded once.
"""
import math
import pathlib
import warnings
from fractions import Fraction

import numpy as np
import scipy
import scipy.linalg as sla
import scipy.signal

OUT = pathlib.Path(__file__).resolve().parent.parent / "tests" / "oracle_cases.hpp"


def rng(seed):
    return np.random.default_rng(seed)


def randn(seed, *shape):
    return rng(seed).standard_normal(shape)


def spd(seed, n):
    M = randn(seed, n, n)
    return M @ M.T / n + np.eye(n)


def graded(n):
    """Diagonal scaling 10^-3 … 10^3: entries of D M D⁻¹ span 10⁻⁶ … 10⁶."""
    return np.diag(np.logspace(-3, 3, n))


def rotation(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s], [s, c]])


def jordan(lam, k):
    return lam * np.eye(k) + np.eye(k, k=1)


def well_conditioned_basis(seed, n):
    """I + 0.3·(random): a non-orthogonal basis with κ(V) of order 1."""
    return np.eye(n) + 0.3 * randn(seed, n, n) / np.sqrt(n)


def zoh(A, B, dt):
    n, m = B.shape
    M = np.zeros((n + m, n + m))
    M[:n, :n] = A * dt
    M[:n, n:] = B * dt
    E = sla.expm(M)
    return E[:n, :n], E[:n, n:]


EPS_FLOAT = float(np.finfo(np.float32).eps)
EPS_DOUBLE = float(np.finfo(np.float64).eps)


def float_policy(c):
    """Run float only while κ·ε_float < 1, i.e. while float can determine a digit of the answer."""
    if c["kappa"] * EPS_FLOAT >= 1.0:
        c["float_ok"] = False
        c["why"] = f"double only: κ·ε_float = {c['kappa'] * EPS_FLOAT:.1g} ≥ 1, float determines no digit of the answer."


def lyap_inv_norm(L):
    """‖L⁻¹‖₂ of an n²×n² Kronecker Lyapunov/Stein operator."""
    return 1.0 / np.linalg.svd(L, compute_uv=False)[-1]


# ---------------------------------------------------------------------------
# expm
# ---------------------------------------------------------------------------

def expm_cases():
    cases = []

    def add(name, A, units=None):
        A = np.asarray(A, dtype=float)
        units = np.ones(A.shape[0]) if units is None else units
        cases.append(dict(name=name, A=A, float_ok=True, why="", units=units))

    add("expm_scalar_1", [[3.7]])
    add("expm_random_3", randn(1, 3, 3) / np.sqrt(3))
    add("expm_random_6", randn(2, 6, 6) / np.sqrt(6))
    add("expm_random_8", 2.0 * randn(3, 8, 8) / np.sqrt(8))
    D = graded(4)
    add("expm_graded_4", D @ (randn(4, 4, 4) / 2.0) @ np.linalg.inv(D), units=np.diag(D))
    add("expm_jordan_4", jordan(-1.0, 4))
    V = well_conditioned_basis(5, 4)
    J = sla.block_diag(jordan(0.3, 3), [[-1.0]])
    add("expm_jordan_similar_4", V @ J @ np.linalg.inv(V))
    V = well_conditioned_basis(6, 5)
    add("expm_repeated_5", V @ np.diag([1.0, 1.0, 1.0, -2.0, -2.0]) @ np.linalg.inv(V))
    S = 50.0 * randn(7, 4, 4) / 2.0
    shift = np.max(np.linalg.eigvals(S).real) - 5.0
    add("expm_large_norm_4", S - shift * np.eye(4))
    V = well_conditioned_basis(8, 3)
    add("expm_stiff_3", V @ np.diag([-1.0, -100.0, -1.0e4]) @ np.linalg.inv(V))
    add("expm_rotation_2", [[0.0, -60.0], [60.0, 0.0]])
    T = np.array([[-1.0, 1.0e4], [0.0, -1.01]])
    add("expm_nonnormal_2", T)
    Rt = rotation(0.3)
    add("expm_nonnormal_rotated_2", Rt.T @ T @ Rt)
    add("expm_nilpotent_5", np.triu(10.0 * randn(9, 5, 5), k=1))
    V = np.eye(3) + 0.9 * np.triu(np.ones((3, 3)), k=1) + 0.05 * randn(10, 3, 3)
    add("expm_near_singular_3", V @ np.diag([1.0e-10, 1.0, -1.0]) @ np.linalg.inv(V))

    for c in cases:
        # Reference and conditioning in natural units M = D⁻¹AD, where the problem is benign.
        D = np.diag(c["units"])
        M = np.linalg.solve(D, c["A"] @ D)
        c["expA"] = D @ sla.expm(M) @ np.linalg.inv(D)
        c["kappa"] = float(sla.expm_cond(M))
        float_policy(c)
    return cases


# ---------------------------------------------------------------------------
# Riccati
# ---------------------------------------------------------------------------

def care_check(A, B, Q, R, N, X):
    K = np.linalg.solve(R, B.T @ X + N.T)
    Acl = A - B @ K
    res = A.T @ X + X @ A - (X @ B + N) @ K + Q
    den = sum(np.linalg.norm(M) for M in (A.T @ X, X @ A, (X @ B + N) @ K, Q))
    margin = float(np.max(np.linalg.eigvals(Acl).real))
    return np.linalg.norm(res) / den, Acl, margin


def dare_check(A, B, Q, R, N, X):
    K = np.linalg.solve(R + B.T @ X @ B, B.T @ X @ A + N.T)
    Acl = A - B @ K
    res = A.T @ X @ A - X - (A.T @ X @ B + N) @ K + Q
    den = sum(np.linalg.norm(M) for M in (A.T @ X @ A, X, (A.T @ X @ B + N) @ K, Q))
    margin = float(np.max(np.abs(np.linalg.eigvals(Acl))))
    return np.linalg.norm(res) / den, Acl, margin


def exact_solve(M, Y):
    """M⁻¹Y in exact rational arithmetic (Gauss–Jordan on object arrays of Fraction)."""
    n = M.shape[0]
    aug = np.concatenate([M, Y], axis=1)
    for c in range(n):
        p = next(r for r in range(c, n) if aug[r, c] != 0)
        aug[[c, p]] = aug[[p, c]]
        aug[c] = aug[c] / aug[c, c]
        for r in range(n):
            if r != c:
                aug[r] = aug[r] - aug[r, c] * aug[c]
    return aug[:, n:]


def res_floor(kind, A, B, Q, R, N, X):
    """Relative residual of the reference X, exactly as stored in double, in units of n·ε.

    The equation is evaluated in exact rational arithmetic on the double data, over the
    same Σ‖terms‖ as care_check / dare_check. It is the residual a correctly rounded answer
    attains: rounding X to double moves the residual by about ε‖Ω‖‖X‖ (Ω the closed-loop
    Lyapunov / Stein operator), which exceeds n·ε·Σ‖terms‖ when ‖A_cl‖ dwarfs the terms
    (cheap control: closed-loop poles near −1/√R).
    """
    A_, B_, Q_, R_, N_, X_ = (np.vectorize(Fraction, otypes=[object])(M) for M in (A, B, Q, R, N, X))
    if kind == "care":
        L = X_ @ B_ + N_
        res = A_.T @ X_ + X_ @ A_ - L @ exact_solve(R_, L.T) + Q_
        den = sum(np.linalg.norm(M) for M in (A.T @ X, X @ A, (X @ B + N) @ np.linalg.solve(R, B.T @ X + N.T), Q))
    else:
        L = A_.T @ X_ @ B_ + N_
        res = A_.T @ X_ @ A_ - X_ - L @ exact_solve(R_ + B_.T @ X_ @ B_, L.T) + Q_
        K = np.linalg.solve(R + B.T @ X @ B, B.T @ X @ A + N.T)
        den = sum(np.linalg.norm(M) for M in (A.T @ X @ A, X, (A.T @ X @ B + N) @ K, Q))
    norm = float(sum(v * v for v in res.ravel())) ** 0.5
    return norm / den / (A.shape[0] * EPS_DOUBLE)


def care_kappa(A, B, Q, R, N, X, Acl):
    """First-order relative condition of X (Byers / Kenney–Hewer form).

    ‖Ω⁻¹‖·(‖Q‖ + 2‖A‖‖X‖ + ‖G‖‖X‖²)/‖X‖ with Ω(Z) = A_clᵀZ + ZA_cl, G = BR⁻¹Bᵀ.
    """
    n = A.shape[0]
    I = np.eye(n)
    L = np.kron(I, Acl.T) + np.kron(Acl.T, I)
    G = B @ np.linalg.solve(R, B.T)
    nX = np.linalg.norm(X)
    num = np.linalg.norm(Q) + 2 * np.linalg.norm(A) * nX + np.linalg.norm(G) * nX**2
    return lyap_inv_norm(L) * num / nX


def dare_kappa(A, B, Q, R, N, X, Acl):
    """First-order relative condition of X: ‖Ω⁻¹‖·(‖Q‖ + 2‖A‖²‖X‖)/‖X‖, Ω(Z) = A_clᵀZA_cl − Z."""
    n = A.shape[0]
    L = np.kron(Acl.T, Acl.T) - np.eye(n * n)
    nX = np.linalg.norm(X)
    num = np.linalg.norm(Q) + 2 * np.linalg.norm(A) ** 2 * nX
    return lyap_inv_norm(L) * num / nX


def riccati_case(name, A, B, Q, R, N=None, units=None):
    A, B, Q, R = (np.atleast_2d(np.asarray(M, dtype=float)) for M in (A, B, Q, R))
    N = np.zeros(B.shape) if N is None else np.atleast_2d(np.asarray(N, dtype=float))
    units = np.ones(A.shape[0]) if units is None else units
    return dict(name=name, A=A, B=B, Q=Q, R=R, N=N, float_ok=True, why="", units=units)


def natural(c):
    """The stored problem in natural units x₀ = D⁻¹x: (D⁻¹AD, D⁻¹B, DQD, R, DN)."""
    D = np.diag(c["units"])
    return np.linalg.solve(D, c["A"] @ D), np.linalg.solve(D, c["B"]), D @ c["Q"] @ D, c["R"], D @ c["N"]


def solve_reference(c, kind):
    A, B, Q, R, N = natural(c)
    if kind == "care":
        X = sla.solve_continuous_are(A, B, Q, R, s=N)
        check, kappa = care_check, care_kappa
    else:
        X = sla.solve_discrete_are(A, B, Q, R, s=N)
        check, kappa = dare_check, dare_kappa
    X = (X + X.T) / 2
    rel, Acl, margin = check(A, B, Q, R, N, X)
    for _ in range(3):
        # Newton (Kleinman / Hewer) refinement: F'(X)[dX] = −F(X), kept while it helps.
        if kind == "care":
            F = A.T @ X + X @ A - (X @ B + N) @ np.linalg.solve(R, B.T @ X + N.T) + Q
            dX = sla.solve_continuous_lyapunov(Acl.T, -F)
        else:
            K = np.linalg.solve(R + B.T @ X @ B, B.T @ X @ A + N.T)
            F = A.T @ X @ A - X - (A.T @ X @ B + N) @ K + Q
            dX = sla.solve_discrete_lyapunov(Acl.T, F)
        Xn = X + (dX + dX.T) / 2
        rel_n, Acl_n, margin_n = check(A, B, Q, R, N, Xn)
        if not rel_n < rel:
            break
        X, rel, Acl, margin = Xn, rel_n, Acl_n, margin_n
    stable = margin < 0 if kind == "care" else margin < 1
    assert stable and rel < 1e-10, f"{c['name']}: scipy reference not trusted (res {rel:g}, margin {margin:g})"
    Di = np.diag(1.0 / c["units"])
    c["X"] = Di @ X @ Di
    c["kappa"] = float(kappa(A, B, Q, R, N, X, Acl))
    c["margin"] = margin
    c["res_floor"] = res_floor(kind, A, B, Q, R, N, X)
    float_policy(c)


def scaled(D, A, B, Q, N):
    """Change of units x = D x₀: A → DAD⁻¹, B → DB, Q → D⁻¹QD⁻¹, N → D⁻¹N."""
    Di = np.linalg.inv(D)
    return D @ A @ Di, D @ B, Di @ Q @ Di, Di @ N


def care_cases():
    ok, fail = [], []
    ok.append(riccati_case("care_scalar_1", [[2.0]], [[1.0]], [[1.0]], [[1.0]]))
    ok.append(riccati_case("care_random_3x1", randn(11, 3, 3), randn(12, 3, 1), np.eye(3), [[1.0]]))
    C = randn(13, 2, 5)
    ok.append(riccati_case("care_random_5x2", randn(14, 5, 5), randn(15, 5, 2), C.T @ C, spd(16, 2)))
    ok.append(riccati_case("care_random_8x3", randn(17, 8, 8) / np.sqrt(8), randn(18, 8, 3), np.eye(8), np.eye(3)))
    A0, B0 = randn(19, 4, 4), randn(20, 4, 1)
    A, B, Q, N = scaled(graded(4), A0, B0, np.eye(4), np.zeros((4, 1)))
    ok.append(riccati_case("care_graded_4", A, B, Q, [[1.0]], N, units=np.diag(graded(4))))
    zeta, w = 1e-4, 10.0
    ok.append(riccati_case("care_lightly_damped_2", [[0.0, 1.0], [-w * w, -2 * zeta * w]], [[0.0], [1.0]], np.eye(2), [[1.0]]))
    e1 = np.eye(4)[:, :1]
    ok.append(riccati_case("care_integrator_chain_4", np.eye(4, k=1), np.eye(4)[:, 3:], e1 @ e1.T, [[1.0]]))
    ok.append(riccati_case("care_near_unstabilizable_3", np.diag([1.0, -1.0, 0.5]), [[1e-3], [1.0], [1.0]], np.eye(3), [[1.0]]))
    ok.append(riccati_case("care_cheap_control_3", randn(21, 3, 3), randn(22, 3, 2), np.eye(3), 1e-8 * np.eye(2)))
    ok.append(riccati_case("care_unstable_4", 3.0 * randn(23, 4, 4) / 2.0 + 4.0 * np.eye(4), randn(24, 4, 2), np.eye(4), np.eye(2)))
    ok.append(riccati_case("care_cross_term_3", randn(25, 3, 3), randn(26, 3, 1), 2.0 * np.eye(3), [[1.0]], 0.5 * randn(27, 3, 1) / np.sqrt(3)))
    ok.append(riccati_case("care_cross_term_4x2", randn(28, 4, 4), randn(29, 4, 2), 3.0 * np.eye(4), spd(30, 2), 0.4 * randn(31, 4, 2)))
    for c in ok:
        solve_reference(c, "care")

    fail.append(riccati_case("care_fail_unstabilizable_2", np.diag([1.0, -1.0]), [[0.0], [1.0]], np.eye(2), [[1.0]]))
    osc = sla.block_diag([[0.0, 1.0], [-1.0, 0.0]], [[-1.0]])
    fail.append(riccati_case("care_fail_imaginary_uncontrollable_3", osc, [[0.0], [0.0], [1.0]], np.eye(3), [[1.0]]))
    return ok, fail


def dare_cases():
    ok, fail = [], []
    ok.append(riccati_case("dare_scalar_1", [[1.5]], [[1.0]], [[1.0]], [[1.0]]))
    A = randn(41, 3, 3)
    ok.append(riccati_case("dare_random_3x1", 1.3 * A / np.max(np.abs(np.linalg.eigvals(A))), randn(42, 3, 1), np.eye(3), [[1.0]]))
    A, C = randn(43, 5, 5), randn(44, 2, 5)
    ok.append(riccati_case("dare_random_5x2", 1.1 * A / np.max(np.abs(np.linalg.eigvals(A))), randn(45, 5, 2), C.T @ C, spd(46, 2)))
    A = randn(47, 8, 8)
    ok.append(riccati_case("dare_random_8x3", 1.2 * A / np.max(np.abs(np.linalg.eigvals(A))), randn(48, 8, 3), np.eye(8), np.eye(3)))
    A0 = randn(49, 4, 4)
    A0 = 1.2 * A0 / np.max(np.abs(np.linalg.eigvals(A0)))
    A, B, Q, N = scaled(graded(4), A0, randn(50, 4, 1), np.eye(4), np.zeros((4, 1)))
    ok.append(riccati_case("dare_graded_4", A, B, Q, [[1.0]], N, units=np.diag(graded(4))))
    zeta, w = 1e-4, 10.0
    Ad, Bd = zoh(np.array([[0.0, 1.0], [-w * w, -2 * zeta * w]]), np.array([[0.0], [1.0]]), 0.01)
    ok.append(riccati_case("dare_lightly_damped_2", Ad, Bd, np.eye(2), [[1.0]]))
    e1 = np.eye(3)[:, :1]
    ok.append(riccati_case("dare_integrator_chain_3", jordan(1.0, 3), np.eye(3)[:, 2:], e1 @ e1.T, [[1.0]]))
    ok.append(riccati_case("dare_near_unstabilizable_3", np.diag([1.2, 0.5, 1.05]), [[1e-3], [1.0], [1.0]], np.eye(3), [[1.0]]))
    ok.append(riccati_case("dare_cheap_control_3", randn(51, 3, 3), randn(52, 3, 2), np.eye(3), 1e-8 * np.eye(2)))
    ok.append(riccati_case("dare_singular_r_2", [[1.0, 0.1], [0.0, 1.0]], [[0.005], [0.1]], np.eye(2), [[0.0]]))
    ok.append(riccati_case("dare_cross_term_3", randn(53, 3, 3) / np.sqrt(3), randn(54, 3, 1), 2.0 * np.eye(3), [[1.0]], 0.5 * randn(55, 3, 1) / np.sqrt(3)))
    ok.append(riccati_case("dare_cross_term_4x2", randn(56, 4, 4) / 2.0, randn(57, 4, 2), 3.0 * np.eye(4), spd(58, 2), 0.4 * randn(59, 4, 2)))
    for c in ok:
        solve_reference(c, "dare")

    fail.append(riccati_case("dare_fail_unstabilizable_2", np.diag([1.5, 0.5]), [[0.0], [1.0]], np.eye(2), [[1.0]]))
    unit = sla.block_diag(rotation(0.3), [[0.5]])
    fail.append(riccati_case("dare_fail_unit_circle_uncontrollable_3", unit, [[0.0], [0.0], [1.0]], np.eye(3), [[1.0]]))
    return ok, fail


# ---------------------------------------------------------------------------
# Shared helpers for the factorization, matrix-function and LTI families
# ---------------------------------------------------------------------------

def make_case(name, **fields):
    return dict(name=name, float_ok=True, why="", **fields)


def frac(M):
    """Exact rational copy of a double array (object dtype)."""
    return np.vectorize(Fraction, otypes=[object])(np.asarray(M, dtype=float))


def to_double(M):
    return np.vectorize(float, otypes=[float])(M)


def exact_det(M):
    """det(M) of a Fraction array by exact Gaussian elimination."""
    M = M.copy()
    n = M.shape[0]
    det = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r, c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            M[[c, p]] = M[[p, c]]
            det = -det
        det *= M[c, c]
        for r in range(c + 1, n):
            M[r] = M[r] - (M[r, c] / M[c, c]) * M[c]
    return det


def exact_solution(A, b):
    """A⁻¹b in exact arithmetic on the double data, rounded once."""
    return to_double(exact_solve(frac(A), frac(np.asarray(b, dtype=float)[:, None])))[:, 0]


def random_orthogonal(seed, n):
    Q, R = np.linalg.qr(randn(seed, n, n))
    return Q * np.sign(np.diag(R))


def similar(seed, J):
    """V·J·V⁻¹ with a well-conditioned basis V, and κ₂(V)."""
    V = well_conditioned_basis(seed, J.shape[0])
    return V @ J @ np.linalg.inv(V), float(np.linalg.cond(V))


def sym(A):
    """Exactly symmetric copy: the lower triangle mirrored."""
    A = np.asarray(A, dtype=float)
    return np.tril(A) + np.tril(A, -1).T


def shifted_stable(seed, n, margin=0.5):
    """A random n×n matrix shifted so its rightmost eigenvalue sits at −margin."""
    A = randn(seed, n, n) / np.sqrt(n)
    return A - (np.max(np.linalg.eigvals(A).real) + margin) * np.eye(n)


def scaled_radius(seed, n, rho):
    """A random n×n matrix scaled to spectral radius rho."""
    A = randn(seed, n, n)
    return rho * A / np.max(np.abs(np.linalg.eigvals(A)))


def oscillator(zeta, w):
    return np.array([[0.0, 1.0], [-w * w, -2.0 * zeta * w]])


def natural_similarity(A, units):
    """D⁻¹AD for D = diag(units)."""
    return A * (1.0 / units)[:, None] * units[None, :]


def ev_cond(A):
    """Eigenvalues of A and their condition numbers 1/|yᴴx| (unit left and right vectors)."""
    w, vl, vr = sla.eig(A, left=True, right=True)
    s = np.abs(np.sum(vl.conj() * vr, axis=0)) / (np.linalg.norm(vl, axis=0) * np.linalg.norm(vr, axis=0))
    return w, 1.0 / s


def cluster_orders(w, clusters):
    """Per eigenvalue: k for the k eigenvalues nearest each (λ, k) Jordan cluster, else 1."""
    order = np.ones(len(w))
    for lam, k in clusters:
        order[np.argsort(np.abs(w - lam))[:k]] = k
    return order


def sorted_roots(w, *per_root):
    """Roots sorted by (real, imag), with arrays that travel with them."""
    idx = np.lexsort((np.round(w.imag, 12), np.round(w.real, 12)))
    return (w[idx],) + tuple(np.asarray(a)[idx] for a in per_root)


def defective_cond(A, lam, k, seed, trials=20):
    """Statistical condition of a k-fold defective eigenvalue: max over random E of (|Δλ|/‖A‖)^k / (‖E‖/‖A‖).

    |Δλ| ≈ (cond·‖E‖/‖A‖)^{1/k}·‖A‖ for a Jordan block of size k; the estimate keeps the
    worst of a few random perturbations of relative size 1e-10.
    """
    nA = np.linalg.norm(A)
    worst = 0.0
    for t in range(trials):
        E = randn(seed + t, *A.shape)
        E *= 1e-10 * nA / np.linalg.norm(E)
        w = np.linalg.eigvals(A + E)
        shift = np.max(np.sort(np.abs(w - lam))[:k])
        worst = max(worst, (shift / nA) ** k / 1e-10)
    return 2.0 * worst


# ---------------------------------------------------------------------------
# LU, QR, Cholesky, LDLᵀ, SVD
# ---------------------------------------------------------------------------

def lu_cases():
    cases = []

    def add(name, A, seed):
        A = np.asarray(A, dtype=float)
        b = randn(seed, A.shape[0])
        sv = np.linalg.svd(A, compute_uv=False)
        cases.append(make_case(name, A=A, b=b, x=exact_solution(A, b), absdet=float(abs(exact_det(frac(A)))),
                          inv_norm=1.0 / sv[-1], kappa=sv[0] / sv[-1]))

    add("lu_random_2", randn(101, 2, 2), 102)
    add("lu_random_4", randn(103, 4, 4), 104)
    add("lu_random_6", randn(105, 6, 6), 106)
    add("lu_needs_pivot_3", [[0.0, 1.0, 2.0], [1.0, 0.0, 3.0], [4.0, -3.0, 8.0]], 107)
    add("lu_hilbert_5", sla.hilbert(5), 108)
    add("lu_graded_rows_4", np.diag(np.logspace(-4, 4, 4)) @ randn(109, 4, 4), 110)
    U, V = random_orthogonal(111, 5), random_orthogonal(112, 5)
    add("lu_near_singular_5", U @ np.diag([1.0, 1.0, 0.5, 0.5, 1e-9]) @ V.T, 113)
    W = np.eye(6) - np.tril(np.ones((6, 6)), -1)
    W[:, -1] = 1.0
    add("lu_wilkinson_growth_6", W, 114)
    for c in cases:
        float_policy(c)
    return cases


def qr_cases():
    cases = []

    def add(name, A):
        A = np.asarray(A, dtype=float)
        _, R = sla.qr(A, mode="economic")
        R = np.where((np.diag(R) < 0)[:, None], -R, R)
        units = np.linalg.norm(A, axis=0)
        cases.append(make_case(name, A=A, R=R, units=units, kappa=float(np.linalg.cond(A / units))))

    add("qr_random_3x3", randn(121, 3, 3))
    add("qr_random_4x2", randn(122, 4, 2))
    add("qr_random_6x4", randn(123, 6, 4))
    add("qr_hilbert_5x5", sla.hilbert(5))
    add("qr_vandermonde_6x4", np.vander(np.linspace(0.0, 1.0, 6), 4, increasing=True))
    M = randn(124, 5, 3)
    M[:, 2] = M[:, 0] + M[:, 1] + 1e-9 * randn(125, 5)
    add("qr_near_rank_deficient_5x3", M)
    add("qr_graded_columns_5x3", randn(126, 5, 3) @ np.diag([1e-4, 1.0, 1e4]))
    for c in cases:
        float_policy(c)
        qr_float_policy(c)
    return cases


def qr_float_policy(c):
    """Run float only while every |R(j, j)| clears qr_decompose's documented rank-loss floor with room.

    qr_decompose() stores R(j, j) = 0 (rank loss) when |R(j, j)| ≤ m·64·ε·‖aⱼ‖ (decomposition.hpp).
    """
    if not c["float_ok"]:
        return
    A, R = c["A"], c["R"]
    rel = np.abs(np.diag(R)) / np.linalg.norm(A, axis=0)
    floor = A.shape[0] * 64 * EPS_FLOAT
    if np.min(rel) <= 4 * floor:
        c["float_ok"] = False
        c["why"] = (f"double only: |R(j, j)| / ‖aⱼ‖ = {np.min(rel):.1g} is within 4× of qr_decompose's "
                    f"rank-loss floor m·64·ε_float = {floor:.1g}.")


def chol_cases():
    ok = []

    def add(name, A, seed):
        A = sym(A)
        units = np.sqrt(np.diag(A))
        b = randn(seed, A.shape[0])
        ok.append(make_case(name, A=A, L=sla.cholesky(A, lower=True), units=units, b=b, x=exact_solution(A, b),
                       kappa=float(np.linalg.cond(A / np.outer(units, units)))))

    add("chol_spd_2", spd(131, 2), 132)
    add("chol_spd_4", spd(133, 4), 134)
    add("chol_spd_6", spd(135, 6), 136)
    add("chol_hilbert_5", sla.hilbert(5), 137)
    add("chol_hilbert_6", sla.hilbert(6), 138)
    add("chol_graded_4", graded(4) @ spd(139, 4) @ graded(4), 140)
    Q = random_orthogonal(141, 4)
    add("chol_near_singular_4", Q @ np.diag([1.0, 1.0, 1.0, 1e-10]) @ Q.T, 142)
    Q = random_orthogonal(143, 5)
    add("chol_repeated_eigenvalues_5", Q @ np.diag([2.0, 2.0, 2.0, 1.0, 1.0]) @ Q.T, 144)
    for c in ok:
        float_policy(c)
        chol_float_policy(c)
    Q = random_orthogonal(145, 3)
    fail = [make_case("chol_fail_indefinite_3", A=sym(Q @ np.diag([1.0, 2.0, -0.5]) @ Q.T))]
    return ok, fail


def chol_float_policy(c):
    """Run float only while every pivot clears cholesky's documented refusal floor with room.

    cholesky() refuses a pivot dᵢ at or below 16·(i+1)·ε times the terms that formed it,
    aᵢᵢ + Σₖ lᵢₖ² (decomposition.hpp); a pivot within 4× of that floor in float tests the
    refusal, not the factor.
    """
    if not c["float_ok"]:
        return
    A, L = c["A"], c["L"]
    for i in range(A.shape[0]):
        pivot = L[i, i] ** 2
        terms = A[i, i] + np.sum(L[i, :i] ** 2)
        floor = 16 * (i + 1) * EPS_FLOAT
        if pivot / terms <= 4 * floor:
            c["float_ok"] = False
            c["why"] = (f"double only: pivot {i} is {pivot / terms:.1g} of its terms, within 4× of "
                        f"cholesky's refusal floor 16·{i + 1}·ε_float = {floor:.1g}.")
            return


def exact_ldlt(A):
    """Unpivoted A = L·diag(d)·Lᵀ of a Fraction array (unique when every leading minor is nonzero)."""
    n = A.shape[0]
    L = np.array([[Fraction(int(i == j)) for j in range(n)] for i in range(n)], dtype=object)
    d = [Fraction(0)] * n
    for j in range(n):
        d[j] = A[j, j] - sum((L[j, k] * L[j, k] * d[k] for k in range(j)), Fraction(0))
        for i in range(j + 1, n):
            L[i, j] = (A[i, j] - sum((L[i, k] * L[j, k] * d[k] for k in range(j)), Fraction(0))) / d[j]
    return L, d


def ldlt_cases():
    ok = []

    def add(name, A, seed, units=None):
        A = sym(A)
        n = A.shape[0]
        units = np.ones(n) if units is None else units
        L, d = exact_ldlt(frac(A))
        A0 = A / np.outer(units, units)
        b = randn(seed, n)
        kappa = max(np.linalg.cond(A0[:k, :k]) for k in range(1, n + 1))
        ok.append(make_case(name, A=A, L=to_double(L), d=np.array([float(v) for v in d]), units=units, b=b,
                       x=exact_solution(A, b), neg=int(np.sum(np.linalg.eigvalsh(A0) < 0)), kappa=float(kappa)))

    add("ldlt_spd_3", spd(151, 3), 152)
    add("ldlt_negative_definite_3", -spd(153, 3), 154)
    Q = random_orthogonal(155, 4)
    add("ldlt_indefinite_4", Q @ np.diag([3.0, -1.0, 2.0, -0.5]) @ Q.T, 156)
    H, Aeq = spd(157, 3), randn(158, 2, 3)
    add("ldlt_kkt_5", np.block([[H, Aeq.T], [Aeq, -0.1 * np.eye(2)]]), 159)
    H, Aeq = spd(160, 4), randn(161, 2, 4)
    add("ldlt_kkt_6", np.block([[H, Aeq.T], [Aeq, np.zeros((2, 2))]]), 162)
    D = np.diag(np.logspace(-2, 2, 5))
    H, Aeq = spd(163, 3), randn(164, 2, 3)
    add("ldlt_graded_kkt_5", D @ np.block([[H, Aeq.T], [Aeq, -0.1 * np.eye(2)]]) @ D, 165, units=np.diag(D))
    for c in ok:
        float_policy(c)
    fail = [
        make_case("ldlt_fail_zero_pivot_2", A=np.array([[0.0, 1.0], [1.0, 0.0]])),
        make_case("ldlt_fail_singular_minor_3", A=np.array([[1.0, 2.0, 0.0], [2.0, 4.0, 1.0], [0.0, 1.0, 1.0]])),
    ]
    return ok, fail


def svd_cases():
    cases = []

    def add(name, A):
        A = np.asarray(A, dtype=float)
        cases.append(make_case(name, A=A, s=sla.svdvals(A)))

    add("svd_random_3x3", randn(171, 3, 3))
    add("svd_random_6x4", randn(172, 6, 4))
    add("svd_random_3x5", randn(173, 3, 5))
    add("svd_hilbert_6", sla.hilbert(6))
    add("svd_rank_2_5x4", randn(174, 5, 2) @ randn(175, 2, 4))
    add("svd_repeated_4", 2.0 * random_orthogonal(176, 4))
    U, V = random_orthogonal(177, 5), random_orthogonal(178, 5)
    add("svd_clustered_5", U @ np.diag([3.0, 1.0 + 1e-10, 1.0, 1.0 - 1e-10, 0.5]) @ V.T)
    add("svd_graded_columns_5x3", randn(179, 5, 3) @ np.diag([1e-4, 1.0, 1e4]))
    return cases


# ---------------------------------------------------------------------------
# Eigenvalues (and the Schur basis, checked on the same cases)
# ---------------------------------------------------------------------------

def eig_cases():
    cases = []

    def add(name, A, units=None, basis_cond=None, clusters=()):
        A = np.asarray(A, dtype=float)
        n = A.shape[0]
        units = np.ones(n) if units is None else units
        w, cond = ev_cond(natural_similarity(A, units))
        if basis_cond is not None:
            # Repeated or defective: per-eigenvalue conditions are undefined; use κ₂(V) of the construction.
            cond = np.full(n, basis_cond)
        w, cond, order = sorted_roots(w, cond, cluster_orders(w, clusters))
        cases.append(make_case(name, A=A, ev_re=w.real, ev_im=w.imag, ev_cond=cond, ev_order=order, units=units,
                          kappa=float(np.max(cond))))

    add("eig_random_2", randn(181, 2, 2))
    add("eig_random_3", randn(182, 3, 3))
    add("eig_random_5", randn(183, 5, 5))
    add("eig_random_6", randn(184, 6, 6))
    add("eig_symmetric_4", sym(randn(185, 4, 4)))
    A, _ = similar(186, sla.block_diag([[-0.5, 2.0], [-2.0, -0.5]], [[0.1, 0.3], [-0.3, 0.1]]))
    add("eig_complex_pairs_4", A)
    A, k = similar(187, np.diag([1.0, 1.0, 1.0, -2.0, -2.0]))
    add("eig_repeated_5", A, basis_cond=k)
    A, k = similar(188, sla.block_diag(jordan(0.5, 3), [[-1.0]]))
    add("eig_jordan_4", A, basis_cond=k, clusters=[(0.5, 3)])
    add("eig_companion_6", sla.companion(np.poly(np.arange(1.0, 7.0))))
    add("eig_graded_4", graded(4) @ (randn(189, 4, 4) / 2.0) @ np.linalg.inv(graded(4)), units=np.diag(graded(4)))
    add("eig_lightly_damped_2", oscillator(1e-4, 10.0))
    R = random_orthogonal(190, 3)
    add("eig_nonnormal_3", R.T @ np.array([[-1.0, 100.0, 0.0], [0.0, -1.01, 100.0], [0.0, 0.0, -1.02]]) @ R)
    for c in cases:
        float_policy(c)
    return cases


# ---------------------------------------------------------------------------
# logm, sqrtm
# ---------------------------------------------------------------------------

def frechet_kron(f, X):
    """n²×n² matrix of the linear map E ↦ f(X, E), column by column."""
    n = X.shape[0]
    K = np.zeros((n * n, n * n))
    for k in range(n * n):
        E = np.zeros(n * n)
        E[k] = 1.0
        K[:, k] = f(X, E.reshape(n, n)).ravel()
    return K


def logm_cases():
    cases = []

    def add(name, A, units=None):
        A = np.asarray(A, dtype=float)
        units = np.ones(A.shape[0]) if units is None else units
        A0 = natural_similarity(A, units)
        X0 = np.real_if_close(sla.logm(A0))
        assert X0.dtype == float and np.linalg.norm(sla.expm(X0) - A0) <= 1e-10 * np.linalg.norm(A0), name
        # L_log(A) inverts L_exp(log A) (Higham 2008, §11.2).
        K = frechet_kron(lambda X, E: sla.expm_frechet(X, E, compute_expm=False), X0)
        kappa = np.linalg.norm(np.linalg.inv(K), 2) * np.linalg.norm(A0) / np.linalg.norm(X0)
        cases.append(make_case(name, A=A, X=natural_similarity(X0, 1.0 / units), units=units, kappa=float(kappa)))

    add("logm_spd_3", spd(201, 3))
    add("logm_exp_random_4", sla.expm(randn(202, 4, 4) / 2.0))
    add("logm_rotation_2", 1.3 * rotation(2.5))
    add("logm_jordan_3", similar(203, jordan(2.0, 3))[0])
    add("logm_near_identity_4", np.eye(4) + 1e-6 * randn(204, 4, 4))
    D = graded(4)
    add("logm_graded_4", D @ sla.expm(randn(205, 4, 4) / 2.0) @ np.linalg.inv(D), units=np.diag(D))
    add("logm_random_6", randn(206, 6, 6) / np.sqrt(6) + 2.0 * np.eye(6))
    Rt = rotation(0.3)
    add("logm_nonnormal_2", Rt.T @ np.array([[1.0, 100.0], [0.0, 2.0]]) @ Rt)
    for c in cases:
        float_policy(c)
    return cases


def sqrtm_cases():
    cases = []

    def add(name, A, units=None):
        A = np.asarray(A, dtype=float)
        n = A.shape[0]
        units = np.ones(n) if units is None else units
        A0 = natural_similarity(A, units)
        X0 = np.real_if_close(sla.sqrtm(A0))
        assert X0.dtype == float and np.linalg.norm(X0 @ X0 - A0) <= 1e-10 * np.linalg.norm(A0), name
        # L_sqrt(A, E) solves X·L + L·X = E (Higham 2008, §6.1).
        I = np.eye(n)
        K = np.kron(I, X0) + np.kron(X0.T, I)
        nA, nX = np.linalg.norm(A0), np.linalg.norm(X0)
        cases.append(make_case(name, A=A, X=natural_similarity(X0, 1.0 / units), units=units,
                          kappa=float(np.linalg.norm(np.linalg.inv(K), 2) * nA / nX), alpha=float(nX * nX / nA)))

    add("sqrtm_spd_4", spd(211, 4))
    add("sqrtm_shifted_random_5", randn(212, 5, 5) / np.sqrt(5) + 2.0 * np.eye(5))
    add("sqrtm_complex_pairs_4", similar(213, sla.block_diag([[1.0, 2.0], [-2.0, 1.0]], [[0.5, 0.2], [-0.2, 0.5]]))[0])
    add("sqrtm_jordan_3", similar(214, jordan(1.5, 3))[0])
    D = graded(4)
    add("sqrtm_graded_4", D @ spd(216, 4) @ np.linalg.inv(D), units=np.diag(D))
    add("sqrtm_large_norm_3", similar(217, np.diag([1e4, 1e5, 1e6]))[0])
    Rt = rotation(0.3)
    add("sqrtm_nonnormal_2", Rt.T @ np.array([[1.0, 100.0], [0.0, 4.0]]) @ Rt)
    for c in cases:
        float_policy(c)
    return cases


# ---------------------------------------------------------------------------
# Lyapunov
# ---------------------------------------------------------------------------

def lyap_case(name, A, Q, discrete, units=None):
    """AX + XAᵀ + Q = 0 or AXAᵀ − X + Q = 0, solved in natural units: A₀ = D⁻¹AD, Q₀ = D⁻¹QD⁻¹, X = DX₀D."""
    A, Q = np.asarray(A, dtype=float), np.asarray(Q, dtype=float)
    n = A.shape[0]
    units = np.ones(n) if units is None else units
    A0 = natural_similarity(A, units)
    Q0 = Q / np.outer(units, units)
    I = np.eye(n)
    if discrete:
        X0 = sla.solve_discrete_lyapunov(A0, Q0)
        L = np.kron(A0, A0) - np.eye(n * n)
        num = 2 * np.linalg.norm(A0) ** 2 * np.linalg.norm(X0) + np.linalg.norm(Q0)
    else:
        X0 = sla.solve_continuous_lyapunov(A0, -Q0)
        L = np.kron(I, A0) + np.kron(A0, I)
        num = 2 * np.linalg.norm(A0) * np.linalg.norm(X0) + np.linalg.norm(Q0)
    c = make_case(name, A=A, Q=Q, X=X0 * np.outer(units, units), units=units,
             kappa=float(lyap_inv_norm(L) * num / np.linalg.norm(X0)))
    float_policy(c)
    lyap_float_policy(c, A0, discrete)
    return c


def lyap_float_policy(c, A0, discrete):
    """Run float only while every Bartels–Stewart pivot clears lyap's documented refusal band with room.

    lyap() / dlyap() refuse a pivot λᵢ + λ̄ⱼ (λᵢλ̄ⱼ − 1) at or below 64·n·ε times smax
    (smax² + 1), smax the largest entry of the Schur form (lyapunov.hpp).
    """
    if not c["float_ok"]:
        return
    S = sla.schur(A0.astype(complex), output="complex")[0]
    lam = np.diag(S)
    smax = np.max(np.abs(S))
    piv = lam[:, None] * lam.conj()[None, :] - 1 if discrete else lam[:, None] + lam.conj()[None, :]
    band = 64 * A0.shape[0] * EPS_FLOAT * (smax * smax + 1 if discrete else smax)
    if np.min(np.abs(piv)) <= 4 * band:
        c["float_ok"] = False
        c["why"] = (f"double only: the smallest Bartels–Stewart pivot {np.min(np.abs(piv)):.1g} is within 4× of "
                    f"the {'dlyap' if discrete else 'lyap'} refusal band {band:.1g} in float.")


def lyap_cases():
    ok = [
        lyap_case("lyap_random_3", shifted_stable(221, 3), np.eye(3), False),
        lyap_case("lyap_random_6", shifted_stable(222, 6), spd(223, 6), False),
        lyap_case("lyap_lightly_damped_2", oscillator(1e-3, 10.0), np.eye(2), False),
        lyap_case("lyap_graded_4", graded(4) @ shifted_stable(224, 4) @ np.linalg.inv(graded(4)),
                  graded(4) @ spd(225, 4) @ graded(4), False, units=np.diag(graded(4))),
        lyap_case("lyap_jordan_4", similar(226, sla.block_diag(jordan(-0.5, 3), [[-2.0]]))[0], np.eye(4), False),
        lyap_case("lyap_nonsymmetric_q_3", shifted_stable(227, 3), randn(228, 3, 3), False),
        lyap_case("lyap_slow_mode_3", similar(229, np.diag([-1e-5, -1.0, -2.0]))[0], np.eye(3), False),
        lyap_case("lyap_gramian_5", shifted_stable(230, 5), randn(231, 5, 1) @ randn(231, 5, 1).T, False),
    ]
    fail = [make_case("lyap_fail_mirrored_2", A=np.diag([1.0, -1.0]), Q=np.eye(2))]
    return ok, fail


def dlyap_cases():
    Ad, _ = zoh(oscillator(1e-3, 10.0), np.array([[0.0], [1.0]]), 0.01)
    ok = [
        lyap_case("dlyap_random_3", scaled_radius(241, 3, 0.8), np.eye(3), True),
        lyap_case("dlyap_random_6", scaled_radius(242, 6, 0.9), spd(243, 6), True),
        lyap_case("dlyap_lightly_damped_2", Ad, np.eye(2), True),
        lyap_case("dlyap_graded_4", graded(4) @ scaled_radius(244, 4, 0.7) @ np.linalg.inv(graded(4)),
                  graded(4) @ spd(245, 4) @ graded(4), True, units=np.diag(graded(4))),
        lyap_case("dlyap_jordan_4", similar(246, sla.block_diag(jordan(0.6, 3), [[-0.3]]))[0], np.eye(4), True),
        lyap_case("dlyap_nonsymmetric_q_3", scaled_radius(247, 3, 0.8), randn(248, 3, 3), True),
        lyap_case("dlyap_unstable_3", similar(249, np.diag([1.5, 0.3, -0.4]))[0], np.eye(3), True),
        lyap_case("dlyap_near_unit_3", similar(250, np.diag([0.99999, 0.5, -0.5]))[0], np.eye(3), True),
    ]
    fail = [make_case("dlyap_fail_reciprocal_2", A=np.diag([2.0, 0.5]), Q=np.eye(2))]
    return ok, fail


# ---------------------------------------------------------------------------
# discretize: ZOH and Tustin
# ---------------------------------------------------------------------------

def zoh_case(name, A, B, Ts, units=None):
    A, B = np.atleast_2d(np.asarray(A, dtype=float)), np.atleast_2d(np.asarray(B, dtype=float))
    n, m = B.shape
    units = np.ones(n) if units is None else units
    A0, B0 = natural_similarity(A, units), B / units[:, None]
    M0 = np.zeros((n + m, n + m))
    M0[:n, :n], M0[:n, n:] = A0 * Ts, B0 * Ts
    Ad0, Bd0 = zoh(A0, B0, Ts)
    c = make_case(name, A=A, B=B, Ts=Ts, Ad=natural_similarity(Ad0, 1.0 / units), Bd=Bd0 * units[:, None], units=units,
             kappa=float(sla.expm_cond(M0)))
    float_policy(c)
    return c


def zoh_cases():
    D = graded(4)
    return [
        zoh_case("zoh_random_3x1", randn(261, 3, 3), randn(262, 3, 1), 0.1),
        zoh_case("zoh_random_6x2", randn(263, 6, 6) / np.sqrt(6), randn(264, 6, 2), 0.05),
        zoh_case("zoh_integrator_chain_4", np.eye(4, k=1), np.eye(4)[:, 3:], 0.5),
        zoh_case("zoh_lightly_damped_2", oscillator(1e-4, 10.0), [[0.0], [1.0]], 0.01),
        zoh_case("zoh_stiff_3", similar(265, np.diag([-1.0, -100.0, -1e4]))[0], randn(266, 3, 1), 1e-3),
        zoh_case("zoh_unstable_2", [[0.0, 1.0], [4.0, 0.0]], [[0.0], [1.0]], 0.1),
        zoh_case("zoh_graded_4", D @ randn(267, 4, 4) @ np.linalg.inv(D), D @ randn(268, 4, 1), 0.1, units=np.diag(D)),
        zoh_case("zoh_long_step_3", randn(269, 3, 3) - 2.0 * np.eye(3), randn(270, 3, 1), 5.0),
    ]


def tustin_case(name, A, B, C, D, Ts, units=None):
    A, B, C, D = (np.atleast_2d(np.asarray(M, dtype=float)) for M in (A, B, C, D))
    n = A.shape[0]
    units = np.ones(n) if units is None else units
    A0, B0, C0 = natural_similarity(A, units), B / units[:, None], C * units[None, :]
    Ad0, Bd0, Cd0, Dd, _ = scipy.signal.cont2discrete((A0, B0, C0, D), Ts, method="bilinear")
    c = make_case(name, A=A, B=B, C=C, D=D, Ts=Ts, Ad=natural_similarity(Ad0, 1.0 / units), Bd=Bd0 * units[:, None],
             Cd=Cd0 / units[None, :], Dd=Dd, units=units, kappa=float(np.linalg.cond(np.eye(n) - A0 * Ts / 2)))
    float_policy(c)
    return c


def tustin_cases():
    G = graded(4)
    ok = [
        tustin_case("tustin_random_3x1", randn(271, 3, 3), randn(272, 3, 1), randn(273, 1, 3), [[0.5]], 0.1),
        tustin_case("tustin_random_5x2", randn(274, 5, 5) / np.sqrt(5), randn(275, 5, 2), randn(276, 2, 5),
                    randn(277, 2, 2), 0.05),
        tustin_case("tustin_lightly_damped_2", oscillator(1e-4, 10.0), [[0.0], [1.0]], [[1.0, 0.0]], [[0.0]], 0.01),
        tustin_case("tustin_stiff_3", similar(278, np.diag([-1.0, -100.0, -1e4]))[0], randn(279, 3, 1),
                    randn(280, 1, 3), [[0.0]], 1e-3),
        tustin_case("tustin_integrator_chain_3", np.eye(3, k=1), np.eye(3)[:, 2:], np.eye(3)[:1], [[0.0]], 0.5),
        tustin_case("tustin_graded_4", G @ randn(281, 4, 4) @ np.linalg.inv(G), G @ randn(282, 4, 1),
                    randn(283, 1, 4) @ np.linalg.inv(G), [[0.0]], 0.1, units=np.diag(G)),
        tustin_case("tustin_near_singular_2", similar(284, np.diag([20.0 * (1.0 - 1e-6), -1.0]))[0], [[1.0], [1.0]],
                    [[1.0, -1.0]], [[0.0]], 0.1),
    ]
    fail = [make_case("tustin_fail_singular_2", A=np.diag([4.0, -1.0]), B=np.array([[1.0], [1.0]]),
                 C=np.array([[1.0, 1.0]]), D=np.array([[0.0]]), Ts=0.5)]
    return ok, fail


# ---------------------------------------------------------------------------
# SISO conversions: ss2tf / ss2zpk, tf2zpk, zpk2tf, tf2ss
# ---------------------------------------------------------------------------

def exact_leverrier(A, B, C, D):
    """num, den of C(sI − A)⁻¹B + D by Faddeev–LeVerrier in exact arithmetic on the double data."""
    n = A.shape[0]
    A_, B_, C_ = frac(A), frac(B), frac(C)
    I = np.array([[Fraction(int(i == j)) for j in range(n)] for i in range(n)], dtype=object)
    M = I * Fraction(0)
    c_prev, den, Ms = Fraction(1), [Fraction(1)], []
    for k in range(1, n + 1):
        M = A_ @ M + I * c_prev
        Ms.append(M)
        c_prev = -np.trace(A_ @ M) / k
        den.append(c_prev)
    num = [Fraction(float(D)) * c for c in den]
    for j in range(n):
        num[j + 1] += (C_ @ Ms[j] @ B_)[0, 0]
    return num, den


def root_scale(coeffs, scale, z, m=1):
    """Size of the ball a root z of multiplicity m moves in, per unit relative coefficient error.

    Σⱼ scaleⱼ·|z|^{deg−j} · m! / |p⁽ᵐ⁾(z)|: a perturbation |Δcⱼ| ≤ η·scaleⱼ moves z by about
    (η·this)^{1/m} (first order for m = 1).
    """
    deg = len(coeffs) - 1
    lift = sum(s * abs(z) ** (deg - j) for j, s in enumerate(scale))
    return float(lift * math.factorial(m) / abs(np.polyval(np.polyder(coeffs, m), z)))


def ss2tf_case(name, A, B, C, D, units=None, basis_cond=None, clusters=()):
    A, B, C = (np.atleast_2d(np.asarray(M, dtype=float)) for M in (A, B, C))
    n = A.shape[0]
    units = np.ones(n) if units is None else units
    A0, B0, C0 = natural_similarity(A, units), B / units[:, None], C * units[None, :]
    num_x, den_x = exact_leverrier(A, B, C, D)
    num, den = np.array([float(v) for v in num_x]), np.array([float(v) for v in den_x])
    sig_ref = scipy.signal.ss2tf(A, B, C, np.array([[D]]))
    assert np.allclose(sig_ref[0][0], num, atol=1e-9 * np.max(np.abs(num))) and np.allclose(sig_ref[1], den), name

    # Coefficient scales: a backward-stable method meets ‖ΔA‖ ≤ η‖A‖, which moves the sum of
    # k×k principal minors by k·C(n,k)·‖A‖ᵏ·η, and C·adj(sI−A)·B by k·C(n,k)·‖A‖ᵏ⁻¹‖B‖‖C‖·η.
    a, bc = np.linalg.norm(A0, 2), np.linalg.norm(B0) * np.linalg.norm(C0)
    k = np.arange(n + 1)
    binom = np.array([float(math.comb(n, j)) for j in k])
    den_scale = k * binom * a**k
    num_scale = k * binom * np.where(k > 0, a ** np.maximum(k - 1, 0), 0.0) * bc + abs(D) * den_scale
    num_scale[0] = abs(D)

    w, cond = ev_cond(A0)
    if basis_cond is not None:
        cond = np.full(n, basis_cond)
    p, p_cond, p_order = sorted_roots(w, cond, cluster_orders(w, clusters))

    if D != 0.0:
        Abar0 = A0 - B0 @ C0 / D
        zw, zcond = ev_cond(Abar0)
        z, zc = sorted_roots(zw, zcond * np.linalg.norm(Abar0))
        gain = float(D)
    else:
        lead = next(j for j in range(1, n + 1) if abs(num[j]) > 1e-8 * num_scale[j]) if np.any(num[1:]) else n + 1
        coeffs = num[lead:]
        zw = np.roots(coeffs) if len(coeffs) > 1 else np.zeros(0, dtype=complex)
        both = np.maximum(num_scale[lead:], np.abs(coeffs))
        z, zc = sorted_roots(zw, [root_scale(coeffs, both, zi) for zi in zw])
        gain = float(num[lead]) if lead <= n else 0.0
    nz = len(z)
    pad = np.zeros(n - nz)
    c = make_case(name, A=A, B=B, C=C, D=float(D), num=num, den=den, num_scale=num_scale, den_scale=den_scale,
             p_re=p.real, p_im=p.imag, p_cond=p_cond, p_order=p_order, z_re=np.concatenate([z.real, pad]),
             z_im=np.concatenate([z.imag, pad]), z_scale=np.concatenate([zc, pad]), n_zeros=nz, gain=gain,
             norm_a=float(np.linalg.norm(A0)), units=units, kappa=float(max(np.max(p_cond), 1.0)))
    float_policy(c)
    return c


def ss2tf_cases():
    G = graded(4)
    Ach = np.eye(4, k=1)
    Ach[3] = [-1.0, -2.0, -3.0, -4.0]
    V = well_conditioned_basis(306, 4)
    Vi = np.linalg.inv(V)
    Aj, kj = similar(305, jordan(-1.0, 3))
    return [
        ss2tf_case("ss2tf_random_2", randn(291, 2, 2), randn(292, 2, 1), randn(293, 1, 2), 0.0),
        ss2tf_case("ss2tf_random_4", randn(294, 4, 4), randn(295, 4, 1), randn(296, 1, 4), 0.0),
        ss2tf_case("ss2tf_random_6", randn(297, 6, 6) / np.sqrt(6), randn(298, 6, 1), randn(299, 1, 6), 0.0),
        ss2tf_case("ss2tf_biproper_3", randn(300, 3, 3), randn(301, 3, 1), randn(302, 1, 3), 0.7),
        ss2tf_case("ss2tf_graded_4", G @ randn(303, 4, 4) @ np.linalg.inv(G), G @ randn(304, 4, 1),
                   randn(304, 1, 4) @ np.linalg.inv(G), 0.0, units=np.diag(G)),
        ss2tf_case("ss2tf_lightly_damped_2", oscillator(1e-4, 10.0), [[0.0], [1.0]], [[1.0, 0.0]], 0.0),
        ss2tf_case("ss2tf_jordan_3", Aj, randn(307, 3, 1), randn(308, 1, 3), 0.0, basis_cond=kj, clusters=[(-1.0, 3)]),
        ss2tf_case("ss2tf_relative_degree_3_4", V @ Ach @ Vi, V @ np.eye(4)[:, 3:], np.array([[1.0, 0.5, 0.0, 0.0]]) @ Vi,
                   0.0),
    ]


def exact_poly(roots, lead=1.0):
    """lead·∏(s − r), highest power first, in exact arithmetic on the double roots (conjugate pairs exact)."""
    c = [(Fraction(lead), Fraction(0))]
    for r in roots:
        rr, ri = Fraction(float(np.real(r))), Fraction(float(np.imag(r)))
        out = []
        for k in range(len(c) + 1):
            a = c[k] if k < len(c) else (Fraction(0), Fraction(0))
            b = c[k - 1] if k >= 1 else (Fraction(0), Fraction(0))
            out.append((a[0] - (rr * b[0] - ri * b[1]), a[1] - (rr * b[1] + ri * b[0])))
        c = out
    assert all(im == 0 for _, im in c)
    return [re for re, _ in c]


def roots_reference(roots, lead):
    """Double coefficients of lead·∏(s − r), and their roots with multiplicities.

    When every coefficient is a double the construction roots are the exact roots; otherwise
    the roots of the rounded polynomial come from numpy (backward stable).
    """
    exact = exact_poly(roots, lead)
    coeffs = np.array([float(v) for v in exact])
    roots = np.asarray(roots, dtype=complex)
    if all(Fraction(float(v)) == v for v in exact):
        ref = roots
    else:
        ref = np.roots(coeffs)
        assert np.allclose(np.sort_complex(ref), np.sort_complex(roots), rtol=1e-6), (roots, ref)
    mult = np.array([float(np.sum(np.abs(roots - r) <= 1e-12 * max(1.0, abs(r)))) for r in ref])
    return coeffs, ref, mult


def tf2zpk_case(name, zeros, poles, gain=1.0, den_lead=1.0):
    num, z, zm = roots_reference(zeros, gain * den_lead)
    den, p, pm = roots_reference(poles, den_lead)

    def scales(coeffs, roots, mult):
        flat = np.full(len(coeffs), np.linalg.norm(coeffs))
        return [root_scale(coeffs, flat, r, int(m)) for r, m in zip(roots, mult)]

    z, zs, zm = sorted_roots(z, scales(num, z, zm), zm)
    p, ps, pm = sorted_roots(p, scales(den, p, pm), pm)
    return make_case(name, num=num, den=den, z_re=z.real, z_im=z.imag, z_scale=zs, z_mult=zm, p_re=p.real, p_im=p.imag,
                p_scale=ps, p_mult=pm, gain=float(gain))


def tf2zpk_cases():
    return [
        tf2zpk_case("tf2zpk_real_3", [-1.0, -3.0], [-2.0, -4.0, -5.0], gain=2.0),
        tf2zpk_case("tf2zpk_complex_pairs_4", [-0.5 + 3j, -0.5 - 3j], [-1 + 2j, -1 - 2j, -0.2 + 2j, -0.2 - 2j], gain=1.5),
        tf2zpk_case("tf2zpk_repeated_pole_3", [-2.0], [-1.0, -1.0, -1.0]),
        tf2zpk_case("tf2zpk_wide_spread_4", [-(2.0**-5)], [-(2.0**-10), -(2.0**-3), -(2.0**3), -(2.0**10)], gain=4.0),
        tf2zpk_case("tf2zpk_biproper_2", [-1 / 3 + 0.4714045207910317j, -1 / 3 - 0.4714045207910317j],
                    [-0.125 + 1.996089927833914j, -0.125 - 1.996089927833914j], gain=3.0),
        tf2zpk_case("tf2zpk_lightly_damped_2", [], [-(2.0**-13) + 1j * np.sqrt(1 - 2.0**-26), -(2.0**-13) - 1j * np.sqrt(1 - 2.0**-26)]),
        tf2zpk_case("tf2zpk_nonmonic_6", [-1.5, -2.5, -3.5], [-1.0, -2.0, -3.0, -4.0, -5.0, -6.0], gain=0.5, den_lead=2.0),
    ]


def zpk2tf_case(name, zeros, poles, gain):
    zeros, poles = np.asarray(zeros, dtype=complex), np.asarray(poles, dtype=complex)

    def esym(roots, lead):
        """Coefficient scales |lead|·eₖ(|r|): the expansion's terms, all of one sign."""
        return np.abs(lead) * np.array([float(v) for v in exact_poly(-np.abs(roots), 1.0)])

    return make_case(name, z_re=zeros.real, z_im=zeros.imag, p_re=poles.real, p_im=poles.imag, gain=float(gain),
                num=np.array([float(v) for v in exact_poly(zeros, gain)]),
                den=np.array([float(v) for v in exact_poly(poles, 1.0)]),
                num_scale=esym(zeros, gain), den_scale=esym(poles, 1.0))


def zpk2tf_cases():
    r = np.abs(randn(311, 5))
    return [
        zpk2tf_case("zpk2tf_real_3", [-1.0, -3.0], [-2.0, -4.0, -5.0], 2.0),
        zpk2tf_case("zpk2tf_complex_4", [-0.5 + 3j, -0.5 - 3j], [-1 + 2j, -1 - 2j, -0.2 + 2j, -0.2 - 2j], 1.5),
        zpk2tf_case("zpk2tf_irrational_5", [-0.3, 0.7 + 1.1j, 0.7 - 1.1j],
                    [-r[0], -r[1] + 1j * r[2], -r[1] - 1j * r[2], -r[3] + 1j * r[4], -r[3] - 1j * r[4]], 0.37),
        zpk2tf_case("zpk2tf_repeated_4", [], [-1.0, -1.0, -1.0, -1.0], 1.0),
        zpk2tf_case("zpk2tf_wide_spread_4", [-1e-3], [-1.1e-3, -0.1, -10.0, -1e3], 7.0),
    ]


def markov(num, den, count):
    """Markov parameters h₀ … h_{count−1} of num/den = Σ hₖ s⁻ᵏ, by exact long division."""
    a = [Fraction(float(v)) for v in den]
    n = len(a) - 1
    b = [Fraction(0)] * (n + 1 - len(num)) + [Fraction(float(v)) for v in num]
    h = []
    for k in range(count):
        bk = b[k] if k <= n else Fraction(0)
        h.append((bk - sum((a[j] * h[k - j] for j in range(1, min(k, n) + 1)), Fraction(0))) / a[0])
    return np.array([float(v) for v in h])


def tf2ss_case(name, num, den):
    num, den = np.asarray(num, dtype=float), np.asarray(den, dtype=float)
    return make_case(name, num=num, den=den, markov=markov(num, den, 2 * (len(den) - 1) + 1))


def tf2ss_cases():
    return [
        tf2ss_case("tf2ss_strictly_proper_3", [1.0, 2.0], [1.0, 6.0, 11.0, 6.0]),
        tf2ss_case("tf2ss_biproper_2", [3.0, 2.0, 1.0], [1.0, 0.25, 4.0]),
        tf2ss_case("tf2ss_nonmonic_4", randn(321, 3), np.concatenate([[2.5], np.abs(randn(322, 4))])),
        tf2ss_case("tf2ss_lightly_damped_2", [100.0], [1.0, 2e-3, 100.0]),
        tf2ss_case("tf2ss_wide_spread_4", [1e3, 1e3], np.poly([-1e-3, -0.1, -10.0, -1e3])),
        tf2ss_case("tf2ss_random_6", randn(323, 4), np.concatenate([[1.7], randn(324, 6)])),
    ]


# ---------------------------------------------------------------------------
# Pole placement
# ---------------------------------------------------------------------------

def place_case(name, A, B, poles, Ts=0.0, seed=0):
    A, B = np.asarray(A, dtype=float), np.atleast_2d(np.asarray(B, dtype=float))
    poles = np.asarray(poles, dtype=complex)
    n, m = B.shape
    distinct = len(np.unique(np.round(poles, 12))) == n
    if distinct or m > 1:
        res = scipy.signal.place_poles(A, B, poles)
        K = res.gain_matrix
        assert np.allclose(np.sort_complex(res.computed_poles), np.sort_complex(poles), atol=1e-8), name
        _, X = np.linalg.eig(A - B @ K)
        cond = np.full(n, float(np.linalg.cond(X / np.linalg.norm(X, axis=0))))
        order = np.ones(n)
    else:
        # A single input with one repeated pole: Ackermann's gain, one Jordan chain.
        lam, k = poles[0].real, n
        A_, b_ = frac(A), frac(B)
        ctrb = np.concatenate([np.linalg.matrix_power(A_, j) @ b_ for j in range(n)], axis=1)
        pA = np.linalg.matrix_power(A_ - np.eye(n, dtype=object) * Fraction(lam), n)
        K = to_double(exact_solve(ctrb.T, np.eye(n, dtype=object)[:, -1:] * Fraction(1)).T @ pA)
        cond = np.full(n, defective_cond(A - B @ K, lam, k, seed))
        order = np.full(n, float(k))
    p, cond, order = sorted_roots(poles, cond, order)
    c = make_case(name, A=A, B=B, Ts=Ts, p_re=p.real, p_im=p.imag, p_cond=cond, p_order=order,
             kappa=float(np.max(cond)))
    float_policy(c)
    return c


def place_cases():
    w = sla.block_diag(oscillator(1e-3, 10.0), oscillator(1e-3, 3.0))
    return [
        place_case("place_siso_3", randn(331, 3, 3), randn(332, 3, 1), [-1.0, -2.0, -3.0]),
        place_case("place_siso_complex_4", randn(333, 4, 4), randn(334, 4, 1), [-1 + 1j, -1 - 1j, -2 + 0.5j, -2 - 0.5j]),
        place_case("place_siso_repeated_3", randn(335, 3, 3), randn(336, 3, 1), [-2.0, -2.0, -2.0], seed=337),
        place_case("place_mimo_4x2", randn(338, 4, 4), randn(339, 4, 2), [-1.0, -2.0, -3.0, -4.0]),
        place_case("place_mimo_complex_6x2", randn(340, 6, 6) / np.sqrt(6), randn(341, 6, 2),
                   [-1 + 2j, -1 - 2j, -2 + 1j, -2 - 1j, -3.0, -4.0]),
        place_case("place_mimo_repeated_4x2", randn(342, 4, 4), randn(343, 4, 2), [-1.0, -1.0, -2.0, -2.0]),
        place_case("place_mimo_5x3", randn(344, 5, 5), randn(345, 5, 3), [-0.5, -1.0, -1.5, -2 + 1j, -2 - 1j]),
        place_case("place_discrete_3x1", scaled_radius(346, 3, 1.2), randn(347, 3, 1), [0.5, 0.6 + 0.2j, 0.6 - 0.2j],
                   Ts=0.1),
        place_case("place_lightly_damped_4x2", w, randn(348, 4, 2), [-5 + 5j, -5 - 5j, -3 + 3j, -3 - 3j]),
    ]


# ---------------------------------------------------------------------------
# Generalized eigenvalues (QZ) and invariant zeros
# ---------------------------------------------------------------------------

def pad_to(v, n):
    v = np.asarray(v)
    return np.concatenate([v, np.zeros(n - len(v), dtype=v.dtype)])


def pencil_eigenvalues(S, N):
    """Generalized eigenvalues of (S, N) (inf where β = 0) and the eigenvectors' factor ‖x‖‖y‖/|yᴴNx|."""
    w, vl, vr = sla.eig(S, N, left=True, right=True, homogeneous_eigvals=True)
    alpha, beta = w
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = np.where(beta != 0, alpha / np.where(beta != 0, beta, 1.0), np.inf)
        vec = np.linalg.norm(vr, axis=0) * np.linalg.norm(vl, axis=0) / np.abs(np.sum(vl.conj() * (N @ vr), axis=0))
    return lam, alpha, beta, vec


def finite_gev(A, B):
    """Finite generalized eigenvalues of (A, B) with conditions (‖A‖ + |λ|‖B‖)·‖x‖‖y‖/|yᴴBx|.

    First order, |Δλ| ≲ η·cond for ‖ΔA‖ ≤ η‖A‖, ‖ΔB‖ ≤ η‖B‖ (Frobenius). An eigenvalue is
    infinite when |β|·‖A‖ ≤ 1e-8·|α|·‖B‖: scipy's β for an exactly singular B is at roundoff.
    """
    lam, alpha, beta, vec = pencil_eigenvalues(A, B)
    nA, nB = np.linalg.norm(A), np.linalg.norm(B)
    keep = np.abs(beta) * nA > 1e-8 * np.abs(alpha) * nB
    return lam[keep], ((nA + np.abs(lam) * nB) * vec)[keep]


def integer_low_rank(seed, n, r):
    """An n×n product of n×r and r×n integer matrices: rank r exactly, every entry a small integer."""
    g = rng(seed)
    M = (g.integers(-3, 4, size=(n, r)) @ g.integers(-3, 4, size=(r, n))).astype(float)
    assert np.linalg.matrix_rank(M) == r
    return M


def qz_case(name, A, B):
    A, B = np.asarray(A, dtype=float), np.asarray(B, dtype=float)
    n = A.shape[0]
    lam, cond = finite_gev(A, B)
    lam, cond = sorted_roots(lam, cond)
    unit = np.linalg.norm(A) / np.linalg.norm(B)
    c = make_case(name, A=A, B=B, ev_re=pad_to(lam.real, n), ev_im=pad_to(lam.imag, n), ev_cond=pad_to(cond, n),
                  n_finite=len(lam), kappa=float(max(np.max(cond, initial=0.0) / unit, 1.0)))
    float_policy(c)
    return c


def qz_cases():
    U, W = well_conditioned_basis(411, 4), well_conditioned_basis(412, 4)
    pairs = sla.block_diag([[-0.5, 2.0], [-2.0, -0.5]], [[0.1, 0.3], [-0.3, 0.1]])
    Q1, Q2 = random_orthogonal(421, 3), random_orthogonal(422, 3)
    H = np.triu(randn(417, 5, 5), -1)
    R = np.triu(randn(418, 5, 5))
    R[0, 0] = 0.0
    R[2, 2] = 0.0
    return [
        qz_case("qz_random_4", randn(401, 4, 4), randn(402, 4, 4)),
        qz_case("qz_identity_b_5", randn(403, 5, 5), np.eye(5)),
        qz_case("qz_spd_b_6", randn(404, 6, 6) / np.sqrt(6), spd(405, 6)),
        qz_case("qz_complex_pairs_4", U @ pairs @ W, U @ W),
        qz_case("qz_singular_b_4", randn(413, 4, 4), integer_low_rank(414, 4, 3)),
        qz_case("qz_rank3_b_5", randn(415, 5, 5), integer_low_rank(416, 5, 3)),
        qz_case("qz_zero_b_diagonal_5", H, R),
        qz_case("qz_all_infinite_3", np.triu(randn(419, 3, 3)) + 3.0 * np.eye(3), np.triu(randn(420, 3, 3), 1)),
        qz_case("qz_nearly_singular_b_3", Q1 @ randn(423, 3, 3) @ Q2, Q1 @ np.diag([1.0, 0.5, 1e-6]) @ Q2),
    ]


def exact_char_poly(M):
    """det(λI − M), highest power first, by Faddeev–LeVerrier on a Fraction array."""
    n = M.shape[0]
    I = np.array([[Fraction(int(i == j)) for j in range(n)] for i in range(n)], dtype=object)
    Mk, c, coeffs = I * Fraction(0), Fraction(1), [Fraction(1)]
    for k in range(1, n + 1):
        Mk = M @ Mk + I * c
        c = -np.trace(M @ Mk) / k
        coeffs.append(c)
    return coeffs


def exact_zero_poly(A, B, C, D):
    """Monic polynomial whose roots are the zeros of a square model with a vector relative degree.

    Falb–Wolovich: rᵢ is the first k at which output i sees an input (Dᵢ ≠ 0: k = 0, else
    Cᵢ·A^(k−1)·B ≠ 0), D*ᵢ that row and C*ᵢ = Cᵢ·A^rᵢ. The decoupled A_z = A − B·D*⁻¹·C* has the
    zeros and Σrᵢ eigenvalues at 0, so det(λI − A_z)/λ^Σrᵢ is the zero polynomial. Exact on the
    stored doubles.
    """
    A_, B_, C_, D_ = frac(A), frac(B), frac(C), frac(D)
    n = A.shape[0]
    Ds, Cs, total = [], [], 0
    for i in range(C.shape[0]):
        row = C_[i:i + 1]
        if any(v != 0 for v in D_[i]):
            Ds.append(D_[i:i + 1])
            Cs.append(row)
            continue
        for r in range(1, n + 1):
            markov = row @ B_
            if any(v != 0 for v in markov.ravel()):
                break
            row = row @ A_
        else:
            raise AssertionError("an output no input reaches")
        Ds.append(markov)
        Cs.append(row @ A_)
        total += r
    Dstar, Cstar = np.concatenate(Ds), np.concatenate(Cs)
    assert exact_det(Dstar) != 0, "no vector relative degree"
    coeffs = exact_char_poly(A_ - B_ @ exact_solve(Dstar.copy(), Cstar))
    assert all(v == 0 for v in coeffs[len(coeffs) - total:])
    return coeffs[:len(coeffs) - total]


def system_pencil(A, B, C, D):
    n = A.shape[0]
    S = np.block([[A, B], [C, D]])
    N = np.zeros_like(S)
    N[:n, :n] = np.eye(n)
    return S, N


def rosenbrock_zeros(A, B, C, D):
    """Finite eigenvalues of a square system pencil, conditions ‖S‖·‖x‖‖y‖/|yᴴNx| (perturbing S only)."""
    S, N = system_pencil(A, B, C, D)
    lam, _, _, vec = pencil_eigenvalues(S, N)
    keep = np.isfinite(lam)
    return lam[keep], (np.linalg.norm(S) * vec)[keep]


def nearest(lam, cond, z, tol):
    i = int(np.argmin(np.abs(lam - z)))
    assert abs(lam[i] - z) <= tol * max(1.0, abs(z)), (z, lam[i])
    return lam[i], cond[i]


def tzero_case(name, A0, B0, C0, D0, zeros=None, ref_D=None, structure=None, alpha=1.0, ux=None, uu=None, uy=None):
    """A zero case in natural units (A0, B0, C0, D0), stored as α·DₓA₀Dₓ⁻¹, α·DₓB₀Dᵤ, D_yC₀Dₓ⁻¹, D_yD₀Dᵤ.

    Square models: the count and identity of the zeros from the exact zero polynomial, the values
    and conditions from scipy's QZ on the system pencil (the root nearest each polynomial root:
    rounding the coefficients would cost more digits than the pencil does). Non-square models:
    the construction zeros; their conditions as zeros of two random square-downs (K·C, K·D or
    B·Kᵀ, D·Kᵀ), whose zeros include them. Each reference zero is certified by
    σ_min(P(z)) ≤ 1e-11·‖S‖. ref_D replaces D for the reference (a D at roundoff level);
    structure = (A, B, C) is the model before a change of coordinates rounded its exact zeros
    (CB = 0) into roundoff-level Markov parameters, for the exact count.
    """
    A0, B0, C0, D0 = (np.atleast_2d(np.asarray(M, dtype=float)) for M in (A0, B0, C0, D0))
    n, m = B0.shape
    p = C0.shape[0]
    Dref = D0 if ref_D is None else np.atleast_2d(np.asarray(ref_D, dtype=float))
    if zeros is None:
        assert m == p
        coeffs = exact_zero_poly(*(structure or (A0, B0, C0)), Dref)
        roots = np.roots([float(v) for v in coeffs]) if len(coeffs) > 1 else np.zeros(0, dtype=complex)
        lam, cond = rosenbrock_zeros(A0, B0, C0, Dref)
        pairs = [nearest(lam, cond, r, 1e-6) for r in roots]
    else:
        pairs = []
        for z in np.asarray(zeros, dtype=complex):
            worst = (z, 0.0)
            for seed in (441, 442):
                if p > m:
                    K = randn(seed, m, p)
                    lam, cond = rosenbrock_zeros(A0, B0, K @ C0, K @ D0)
                else:
                    K = randn(seed, p, m)
                    lam, cond = rosenbrock_zeros(A0, B0 @ K.T, C0, D0 @ K.T)
                _, c = nearest(lam, cond, z, 1e-8)
                worst = (z, max(worst[1], c))
            pairs.append(worst)
    z = np.array([v for v, _ in pairs], dtype=complex)
    zc = np.array([c for _, c in pairs], dtype=float)
    S0, N0 = system_pencil(A0, B0, C0, D0)
    for zi in z:
        smin = np.linalg.svd(S0 - zi * N0, compute_uv=False)[-1]
        assert smin <= 1e-11 * np.linalg.norm(S0), (name, zi, smin)
    z, zc = sorted_roots(z, zc)
    ux, uu, uy = (np.ones(k) if u is None else np.asarray(u, dtype=float) for u, k in ((ux, n), (uu, m), (uy, p)))
    Dx, Du, Dy = np.diag(ux), np.diag(uu), np.diag(uy)
    A = alpha * Dx @ A0 @ np.diag(1.0 / ux)
    B = alpha * Dx @ B0 @ Du
    C = Dy @ C0 @ np.diag(1.0 / ux)
    D = Dy @ D0 @ Du
    scale = np.maximum(np.abs(z), np.linalg.norm(A0))
    c = make_case(name, A=A, B=B, C=C, D=D, z_re=pad_to(alpha * z.real, n), z_im=pad_to(alpha * z.imag, n),
                  z_scale=pad_to(alpha * zc, n), n_zeros=len(z), alpha=float(alpha), units_x=ux, units_u=uu, units_y=uy,
                  kappa=float(max(np.max(zc / scale, initial=0.0), 1.0)))
    float_policy(c)
    return c


def series(first, second):
    """second ∘ first: the outputs of first drive the inputs of second (both (A, B, C, D))."""
    A1, B1, C1, D1 = first
    A2, B2, C2, D2 = second
    n1, n2 = A1.shape[0], A2.shape[0]
    A = np.block([[A1, np.zeros((n1, n2))], [B2 @ C1, A2]])
    return A, np.vstack([B1, B2 @ D1]), np.hstack([D2 @ C1, C2]), D2 @ D1


def hidden(seed, sys):
    """The model in the coordinates of a random orthogonal Q (z = Qᵀx): the structure is not visible."""
    A, B, C, D = sys
    Q = random_orthogonal(seed, A.shape[0])
    return Q.T @ A @ Q, Q.T @ B, C @ Q, D


def first_order_zero(z0, a, width=1, channel=0):
    """diag((s − z₀)/(s + a) on one channel, 1 on the others): one state, D = I."""
    B = np.zeros((1, width))
    B[0, channel] = 1.0
    C = np.zeros((width, 1))
    C[channel, 0] = -(z0 + a)
    return np.array([[-a]]), B, C, np.eye(width)


def nmp_square(seed0):
    """A random 5-state, 2×2 model with a zero in the right half plane (first seed that has one)."""
    for s in range(seed0, seed0 + 50):
        A, B, C = randn(s, 5, 5) / np.sqrt(5), randn(s + 1000, 5, 2), randn(s + 2000, 2, 5)
        lam, _ = rosenbrock_zeros(A, B, C, np.zeros((2, 2)))
        if np.any(lam[np.abs(lam) < 1e3].real > 0.1):
            return A, B, C
    raise AssertionError("no non-minimum-phase draw")


def tzero_cases():
    A_siso, B_siso, C_siso, D_siso = scipy.signal.tf2ss(np.poly([2.0, -3.0]), np.poly([-1.0, -2.0, -4.0, -0.5]))
    V = well_conditioned_basis(431, 4)
    Vi = np.linalg.inv(V)
    siso_nmp = (V @ A_siso @ Vi, V @ B_siso, C_siso @ Vi, D_siso)
    A3, B3, C3, _ = scipy.signal.tf2ss([1.0, 1.0], np.poly([-1.0, -2.0, -3.0]))
    V3 = well_conditioned_basis(432, 3)
    V3i = np.linalg.inv(V3)
    small_d = (V3 @ A3 @ V3i, V3 @ B3, C3 @ V3i)
    A4, B4, C4, _ = scipy.signal.tf2ss([1.0, 2.0], np.poly([-1.0, -2.0, -3.0, -5.0]))
    sq = (shifted_stable(433, 4), randn(434, 4, 2), randn(435, 2, 4), np.zeros((2, 2)))
    mixed_A = randn(436, 5, 5) / np.sqrt(5)
    mixed_B = np.vstack([randn(437, 2, 2), np.zeros((3, 2))])
    mixed_C = np.vstack([randn(438, 1, 5), np.hstack([np.zeros((1, 2)), randn(439, 1, 3)])])
    chains_A = np.array([[0.0, 1.0, 0.0, 0.0], [-2.0, -0.5, 1.0, 0.3], [0.0, 0.0, 0.0, 1.0], [0.4, -1.0, -3.0, -0.8]])
    chains_B = np.array([[0.0, 0.0], [1.0, 0.5], [0.0, 0.0], [0.3, 1.0]])
    chains_C = np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]])
    h_tall = (shifted_stable(445, 2), randn(446, 2, 1), randn(447, 2, 2), np.zeros((2, 1)))
    h_tall_d = (shifted_stable(448, 2), randn(449, 2, 1), randn(450, 3, 2), randn(451, 3, 1))
    h_wide = (shifted_stable(452, 3), randn(453, 3, 3), randn(454, 2, 3), np.zeros((2, 3)))
    unobs_A = sla.block_diag([[-1.5]], shifted_stable(455, 3))
    unobs_A[0, 1:] = randn(456, 3)
    unobs = (unobs_A, randn(457, 4, 2), np.hstack([np.zeros((3, 1)), randn(458, 3, 3)]), np.zeros((3, 2)))
    return [
        tzero_case("tzero_siso_nmp_4", *siso_nmp, structure=(A_siso, B_siso, C_siso)),
        tzero_case("tzero_siso_relative_degree_3_4", A4, B4, C4, [[0.0]]),
        tzero_case("tzero_siso_small_d_3", *small_d, [[1e-6]], structure=(A3, B3, C3)),
        tzero_case("tzero_siso_roundoff_d_3", *small_d, [[1e-20]], ref_D=[[0.0]], structure=(A3, B3, C3)),
        tzero_case("tzero_square_2x2_4", *sq),
        tzero_case("tzero_square_biproper_3x3_5", shifted_stable(424, 5), randn(425, 5, 3), randn(426, 3, 5),
                   randn(427, 3, 3)),
        tzero_case("tzero_square_nmp_2x2_5", *nmp_square(460), np.zeros((2, 2))),
        tzero_case("tzero_mixed_relative_degree_2x2_5", mixed_A, mixed_B, mixed_C, np.zeros((2, 2))),
        tzero_case("tzero_no_finite_zeros_2x2_4", chains_A, chains_B, chains_C, np.zeros((2, 2))),
        tzero_case("tzero_units_2x2_4", *sq, alpha=1e3, ux=np.logspace(-2, 2, 4), uu=[1e4, 1e-3], uy=[1e-3, 1e2]),
        tzero_case("tzero_tall_2x1_3", *hidden(461, series(first_order_zero(2.0, 1.5), h_tall)), zeros=[2.0]),
        tzero_case("tzero_tall_feedthrough_3x1_3", *hidden(462, series(first_order_zero(-0.7, 2.0), h_tall_d)),
                   zeros=[-0.7]),
        tzero_case("tzero_wide_2x3_4", *hidden(463, series(h_wide, first_order_zero(-0.5, 3.0, 2, 0))), zeros=[-0.5]),
        tzero_case("tzero_tall_unobservable_3x2_4", *hidden(464, unobs), zeros=[-1.5]),
    ]


# ---------------------------------------------------------------------------
# Emit
# ---------------------------------------------------------------------------

def lit(x):
    s = repr(float(x))
    assert np.isfinite(x), s
    return s


def arr(name, M):
    vals = ", ".join(lit(v) for v in np.asarray(M, dtype=float).ravel())
    return f"    static constexpr std::array<double, {np.asarray(M).size}> {name}{{{vals}}};"


def emit_case(c, fields, sizes, scalars):
    lines = []
    if c["why"]:
        lines.append(f"// {c['why']}")
    lines.append(f"struct {c['name']} {{")
    lines.append(f'    static constexpr const char* name = "{c["name"]}";')
    for k, v in sizes.items():
        lines.append(f"    static constexpr std::size_t {k} = {v};")
    lines.append(f"    static constexpr bool float_ok = {'true' if c['float_ok'] else 'false'};")
    for k in scalars:
        lines.append(f"    static constexpr double {k} = {lit(c[k])};")
    for k in fields:
        lines.append(arr(k, c[k]))
    lines.append("};")
    return "\n".join(lines)


def emit_tuple(alias, cases):
    names = ", ".join(c["name"] for c in cases)
    return f"using {alias} = std::tuple<{names}>;"


def sq(key):
    return lambda c: {"n": c[key].shape[0]}


def rect(c):
    return {"rows": c["A"].shape[0], "cols": c["A"].shape[1]}


def ss_sizes(c):
    return {"nx": c["B"].shape[0], "nu": c["B"].shape[1]}


def tustin_sizes(c):
    return {"nx": c["B"].shape[0], "nu": c["B"].shape[1], "ny": c["C"].shape[0]}


def tf_sizes(c):
    return {"nnum": len(c["num"]), "nden": len(c["den"])}


def zpk_sizes(c):
    return {"nz": len(c["z_re"]), "np": len(c["p_re"])}


def families():
    """(tuple alias, cases, array fields, sizes, double scalars) for every family after the Riccati ones."""
    chol_ok, chol_fail = chol_cases()
    ldlt_ok, ldlt_fail = ldlt_cases()
    lyap_ok, lyap_fail = lyap_cases()
    dlyap_ok, dlyap_fail = dlyap_cases()
    tustin_ok, tustin_fail = tustin_cases()
    roots = ["p_re", "p_im", "p_cond", "p_order"]
    return [
        ("lu_cases", lu_cases(), ["A", "b", "x"], sq("A"), ["kappa", "inv_norm", "absdet"]),
        ("qr_cases", qr_cases(), ["A", "R", "units"], rect, ["kappa"]),
        ("chol_cases", chol_ok, ["A", "L", "units", "b", "x"], sq("A"), ["kappa"]),
        ("chol_fail_cases", chol_fail, ["A"], sq("A"), []),
        ("ldlt_cases", ldlt_ok, ["A", "L", "d", "units", "b", "x"], lambda c: {"n": c["A"].shape[0], "neg": c["neg"]},
         ["kappa"]),
        ("ldlt_fail_cases", ldlt_fail, ["A"], sq("A"), []),
        ("svd_cases", svd_cases(), ["A", "s"], rect, []),
        ("eig_cases", eig_cases(), ["A", "ev_re", "ev_im", "ev_cond", "ev_order", "units"], sq("A"), ["kappa"]),
        ("logm_cases", logm_cases(), ["A", "X", "units"], sq("A"), ["kappa"]),
        ("sqrtm_cases", sqrtm_cases(), ["A", "X", "units"], sq("A"), ["kappa", "alpha"]),
        ("lyap_cases", lyap_ok, ["A", "Q", "X", "units"], sq("A"), ["kappa"]),
        ("lyap_fail_cases", lyap_fail, ["A", "Q"], sq("A"), []),
        ("dlyap_cases", dlyap_ok, ["A", "Q", "X", "units"], sq("A"), ["kappa"]),
        ("dlyap_fail_cases", dlyap_fail, ["A", "Q"], sq("A"), []),
        ("zoh_cases", zoh_cases(), ["A", "B", "Ad", "Bd", "units"], ss_sizes, ["Ts", "kappa"]),
        ("tustin_cases", tustin_ok, ["A", "B", "C", "D", "Ad", "Bd", "Cd", "Dd", "units"], tustin_sizes, ["Ts", "kappa"]),
        ("tustin_fail_cases", tustin_fail, ["A", "B", "C", "D"], tustin_sizes, ["Ts"]),
        ("ss2tf_cases", ss2tf_cases(),
         ["A", "B", "C", "num", "den", "num_scale", "den_scale"] + roots + ["z_re", "z_im", "z_scale", "units"],
         lambda c: {"nx": c["A"].shape[0], "n_zeros": c["n_zeros"]}, ["D", "gain", "norm_a", "kappa"]),
        ("tf2zpk_cases", tf2zpk_cases(),
         ["num", "den", "z_re", "z_im", "z_scale", "z_mult", "p_re", "p_im", "p_scale", "p_mult"], tf_sizes, ["gain"]),
        ("zpk2tf_cases", zpk2tf_cases(), ["z_re", "z_im", "p_re", "p_im", "num", "den", "num_scale", "den_scale"],
         zpk_sizes, ["gain"]),
        ("tf2ss_cases", tf2ss_cases(), ["num", "den", "markov"], tf_sizes, []),
        ("place_cases", place_cases(), ["A", "B"] + roots, ss_sizes, ["Ts", "kappa"]),
        ("qz_cases", qz_cases(), ["A", "B", "ev_re", "ev_im", "ev_cond"],
         lambda c: {"n": c["A"].shape[0], "n_finite": c["n_finite"]}, ["kappa"]),
        ("tzero_cases", tzero_cases(), ["A", "B", "C", "D", "z_re", "z_im", "z_scale", "units_x", "units_u", "units_y"],
         lambda c: {"nx": c["B"].shape[0], "nu": c["B"].shape[1], "ny": c["C"].shape[0], "n_zeros": c["n_zeros"]},
         ["alpha", "kappa"]),
    ]


def main():
    # Newton refinement steps may be ill-conditioned; each is kept only if the residual drops.
    warnings.simplefilter("ignore", sla.LinAlgWarning)
    ex = expm_cases()
    care_ok, care_fail = care_cases()
    dare_ok, dare_fail = dare_cases()
    more = families()

    def rsizes(c):
        return {"nx": c["B"].shape[0], "nu": c["B"].shape[1]}

    riccati_in = ["A", "B", "Q", "R", "N"]
    parts = [
        "// Copyright Paul Guénette and Damp contributors.",
        "// SPDX-License-Identifier: BSL-1.0",
        "// Distributed under the Boost Software License, Version 1.0.",
        "// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)",
        "",
        "// Generated by tools/gen_oracle.py — do not edit by hand. Re-run the script instead.",
        f"// numpy {np.__version__}, scipy {scipy.__version__}. Matrices are row-major.",
        "//",
        "// expm cases:  A, expA = D·scipy.linalg.expm(D⁻¹AD)·D⁻¹, kappa = scipy.linalg.expm_cond(D⁻¹AD)",
        "//              (relative Frobenius condition number of exp, in natural units).",
        "// care / dare: A, B, Q, R, N (cross term), X = scipy.linalg.solve_*_are(..., s=N) in natural",
        "//              units, then up to 3 Newton steps kept while the residual drops,",
        "//              kappa = first-order relative condition estimate of X,",
        "//              margin = max Re λ (care) or max |λ| (dare) of A − BK at the scipy X,",
        "//              res_floor = relative residual of the scipy X in exact arithmetic, in",
        "//              units of n·ε: what a correctly rounded answer attains.",
        "// units = diag(D): the natural-units state is x₀ = D⁻¹x (all ones unless graded). kappa,",
        "//              margin and the test checks are taken in natural units, so a case that",
        "//              only mixes units is judged like the benign problem underneath.",
        "// *_fail_cases have no stabilizing solution: the design must report failure.",
        "// float_ok = false marks a case run in double only (reason above the struct).",
        "//",
        "// The families below ('exact' = rational arithmetic on the stored doubles, rounded once):",
        "// lu:     b, x = exact A⁻¹b, absdet = exact |det A|, inv_norm = ‖A⁻¹‖₂, kappa = κ₂(A).",
        "// qr:     R = scipy R with its rows signed so diag(R) ≥ 0 (unique for full column rank),",
        "//         units = column norms of A, kappa = κ₂ of A with unit columns.",
        "// chol:   L = scipy cholesky (lower), units = √diag(A), kappa = κ₂(D⁻¹AD⁻¹), x = exact A⁻¹b.",
        "// ldlt:   L, d = exact unpivoted LDLᵀ, neg = negative eigenvalues of A (inertia),",
        "//         kappa = max κ₂ of the leading principal blocks (natural units), x = exact A⁻¹b.",
        "// svd:    s = scipy.linalg.svdvals(A).",
        "// eig:    eigenvalues sorted by (re, im) with per-eigenvalue condition 1/|yᴴx| in natural",
        "//         units, or κ₂ of the construction basis for repeated and defective spectra;",
        "//         ev_order = size of the Jordan block the eigenvalue belongs to.",
        "// logm / sqrtm: X = scipy logm / sqrtm in natural units, kappa = relative Frobenius",
        "//         condition from the Kronecker form of the Fréchet derivative, alpha = ‖X‖²/‖A‖.",
        "// lyap / dlyap: X = scipy solve_continuous_lyapunov(A, −Q) / solve_discrete_lyapunov(A, Q)",
        "//         in natural units (A₀ = D⁻¹AD, Q₀ = D⁻¹QD⁻¹, X = DX₀D), kappa = ‖Ω⁻¹‖ times the",
        "//         size of the equation's terms over ‖X‖ (Ω the Lyapunov / Stein operator).",
        "// zoh:    Ad, Bd from expm([A B; 0 0]·Ts), kappa = expm_cond of that block (natural units).",
        "// tustin: scipy.signal.cont2discrete(..., 'bilinear'), kappa = κ₂(I − A·Ts/2).",
        "// ss2tf:  num, den = exact Faddeev–LeVerrier; num_scale, den_scale = what a backward-stable",
        "//         method can move each coefficient by, per unit of relative error; poles and",
        "//         zeros as eig (zeros: roots of the numerator, or eig(A − BC/D) when D ≠ 0) with",
        "//         z_scale = the zero's sensitivity to that coefficient error.",
        "// tf2zpk: roots of num, den (exact when the coefficients are doubles), multiplicity and",
        "//         sensitivity to a normwise coefficient error ‖Δc‖ ≤ η‖c‖.",
        "// zpk2tf: num, den = exact expansion; scales = |k|·eₖ(|roots|), the size of the terms.",
        "// tf2ss:  markov = h₀ … h₂ₙ of num/den by exact long division (realization-independent).",
        "// place:  requested poles with κ₂ of scipy's closed-loop eigenvectors (Bauer–Fike), or for",
        "//         a single input with one repeated pole a statistical estimate of the Jordan",
        "//         block's sensitivity (p_order = block size).",
        "// qz:     finite generalized eigenvalues of (A, B) from scipy.linalg.eig(A, B) (infinite:",
        "//         |β|‖A‖ ≤ 1e-8|α|‖B‖), ev_cond = (‖A‖ + |λ|‖B‖)‖x‖‖y‖/|yᴴBx|, the first-order",
        "//         sensitivity to a normwise relative perturbation of A and B.",
        "// tzero:  invariant zeros. Square models: count from the exact zero polynomial",
        "//         det(λI − A_z)/λ^Σrᵢ (Falb–Wolovich decoupling, rational arithmetic), values and",
        "//         z_scale = ‖S‖‖x‖‖y‖/|yᴴNx| from scipy's QZ on the system pencil S − λN. Non-square",
        "//         models: the zeros of the construction (a common zero of every output, an",
        "//         unobservable mode) behind a random orthogonal change of state coordinates,",
        "//         z_scale the larger of their conditions in two random square-downs. Every",
        "//         reference zero has σ_min(S − zN) ≤ 1e-11‖S‖. Stored data are α·DₓA₀Dₓ⁻¹,",
        "//         α·DₓB₀Dᵤ, D_yC₀Dₓ⁻¹, D_yD₀Dᵤ for the natural (A₀, B₀, C₀, D₀) (units = diag D,",
        "//         α a time scale): the zeros are α·z₀.",
        "#pragma once",
        "",
        "#include <array>",
        "#include <cstddef>",
        "#include <tuple>",
        "",
        "// clang-format off",
        "namespace oracle {",
        "",
    ]
    for c in ex:
        parts.append(emit_case(c, ["A", "expA", "units"], {"n": c["A"].shape[0]}, ["kappa"]))
        parts.append("")
    for c in care_ok + dare_ok:
        parts.append(emit_case(c, riccati_in + ["X", "units"], rsizes(c), ["kappa", "margin", "res_floor"]))
        parts.append("")
    for c in care_fail + dare_fail:
        parts.append(emit_case(c, riccati_in + ["units"], rsizes(c), []))
        parts.append("")
    for _, cases, fields, sizes, scalars in more:
        for c in cases:
            parts.append(emit_case(c, fields, sizes(c), scalars))
            parts.append("")
    parts += [
        emit_tuple("expm_cases", ex),
        emit_tuple("care_cases", care_ok),
        emit_tuple("care_fail_cases", care_fail),
        emit_tuple("dare_cases", dare_ok),
        emit_tuple("dare_fail_cases", dare_fail),
    ]
    parts += [emit_tuple(alias, cases) for alias, cases, *_ in more]
    parts += [
        "",
        "} // namespace oracle",
        "// clang-format on",
        "",
    ]
    OUT.write_text("\n".join(parts), encoding="utf-8", newline="\n")
    counts = ", ".join(f"{len(cases)} {alias[:-6]}" for alias, cases, *_ in more)
    print(f"wrote {OUT} ({len(ex)} expm, {len(care_ok)}+{len(care_fail)} care, {len(dare_ok)}+{len(dare_fail)} dare, {counts})")


if __name__ == "__main__":
    main()
