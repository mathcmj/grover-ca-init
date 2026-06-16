"""
Two-panel figure for resource analysis

Panel (a): n=100, mu=10, vary k (# constraint sets), |S|=1
  - Cardinality (Prop. 2): Dicke |D^mu_1> per register
  - Parity     (Prop. 4): GHZ_mu_^(X) per register

Panel (b): n=100, k=1, vary mu = 2..20, |S|=1
  - Cardinality (Prop. 3): Dicke |D^mu_1>
  - Parity     (Prop. 4): GHZ_mu_^(X)

y-axis: S_sigma + (D + 2*S_sigma) * kappa_opt  [total gate count, oracle excluded]
ancilla-free MCX decomposition with {'rz', 'h', 'cx'}.
"""

import numpy as np
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import Decompose

# MCX is the only multi-qubit gate we decompose (ancilla-free → one & two-qubit gates).
# All other gates (CRY, CX, H, RZ, X, Z, …) are already ≤2-qubit and count as 1.
MCX_BASIS = ['rz', 'h', 'cx']   # target basis for MCX decomposition

def mcx_noancilla_as_instrs(num_ctrls: int, basis=MCX_BASIS):

    if num_ctrls < 1:
        raise ValueError("num_ctrls must be >= 1")
    k = num_ctrls

    sub = QuantumCircuit(k + 1)
    sub.mcx(list(range(k)), k)

    sub_u = transpile(sub, basis_gates=basis, optimization_level=0)

    while True:
        multi = {ci.operation.name for ci in sub_u.data if ci.operation.num_qubits > 2}
        if not multi:
            break
        pm = PassManager([Decompose(g) for g in sorted(multi)])
        sub_u = pm.run(sub_u)

    assert all(ci.operation.num_qubits <= 2 for ci in sub_u.data), sub_u

    q2i = {q: i for i, q in enumerate(sub_u.qubits)}

    instrs = []
    for ci in sub_u.data:
        op = ci.operation
        qidx = [q2i[q] for q in ci.qubits]
        instrs.append((op, qidx))
    return instrs

def transpiled_mcx(qc: QuantumCircuit, controls: list[int], target: int):
    k = len(controls)
    instrs = mcx_noancilla_as_instrs(k)
    mapping = {i: qc.qubits[q] for i, q in enumerate(controls + [target])}

    for op, local_qidx in instrs:
        mapped_qargs = [mapping[i] for i in local_qidx]
        qc.append(op, mapped_qargs)


# ── circuit builders ─────────────────────────────────────────────────────────

def build_D(n: int) -> QuantumCircuit:
    """2|0><0|^n - I, ancilla-free MCX."""
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.x(i)
    qc.h(n-1)
    transpiled_mcx(qc, list(range(n-1)), n-1)
    qc.h(n-1)
    for i in range(1, n):
        qc.x(i)
    qc.z(0)
    qc.x(0)
    qc.z(0)
    return qc


def scs_n_1(circ, m, i):
    circ.cry(2*np.arccos(np.sqrt(1/m)), control_qubit=i+1, target_qubit=i)
    circ.cx(i, i+1)

def unitary_W(circ, q_list):
    n = len(q_list)
    circ.x(q_list[-1])
    for it in range(n - 1):
        scs_n_1(circ, n-it, q_list[n-2-it])

def build_dicke(mu: int) -> QuantumCircuit:
    qc = QuantumCircuit(mu)
    unitary_W(qc, list(range(mu)))
    return qc


def build_ghz(mu: int) -> QuantumCircuit:
    qc = QuantumCircuit(mu)
    qc.h(0)
    for i in range(mu-1):
        qc.cx(i, i+1)
    qc.z(0)
    for i in range(mu):
        qc.h(i)
    return qc


D_100 = len(build_D(100).data)
print(f"D(100) total gates: {len(build_D(100).data)}")
print(f"D(100) depth:       {build_D(100).depth()}")
print(f"D(100) gate breakdown: {dict(build_D(100).count_ops())}")

D_10 = len(build_D(10).data)
print(f"D(10) total gates: {len(build_D(10).data)}")
print(f"D(10) depth:       {build_D(10).depth()}")
print(f"D(10) gate breakdown: {dict(build_D(10).count_ops())}")


# ── Dicke & GHZ gate counts ──────────────────────────────────────────────────

print("=== Dicke state |D^mu_1> gate counts ===")
dicke_costs = {}
for mu in range(2, 21):
    c = len(build_dicke(mu).data)   # CRY, CX, X each count as 1
    dicke_costs[mu] = c
print("  " + "  ".join(f"mu={mu}: {dicke_costs[mu]}" for mu in [2,5,10,15,20]))

print("\n=== GHZ-type state gate counts ===")
ghz_costs = {}
for mu in range(2, 21):
    c = len(build_ghz(mu).data)     # H, CX each count as 1
    ghz_costs[mu] = c
print("  " + "  ".join(f"mu={mu}: {ghz_costs[mu]}" for mu in [2,5,10,15,20]))


# ── resource formula ─────────────────────────────────────────────────────────

def kappa_opt(F_size: float, S_size: int = 1) -> float:
    if F_size <= S_size:
        return 0.0
    x = np.sqrt(S_size/F_size)
    return max(0.0, round(np.pi/(4.0*np.arcsin(x)) - 0.5))

def R(S_cost: float, D_cost: float, kappa: float) -> float:
    """S + (D + 2S) * kappa  (oracle excluded)."""
    return S_cost + (D_cost + 2.0*S_cost) * kappa


# ── Panel (a): n=100, mu=10, vary k ─────────────────────────────────────────

n = 100
mu_a = 10
k_vals = list(range(11))
R_card, R_par = [], []
for k in k_vals:
    # Cardinality  |F| = mu^k * 2^(n-mu*k)
    F_c = float(mu_a**k) * float(2**(n-mu_a*k))
    S_c = k*dicke_costs[mu_a] + (n-mu_a*k)
    R_card.append(R(S_c, D_100, kappa_opt(F_c)))

    # Parity  |F| = 2^(n-k)  (one parity eq. halves search space)
    F_p = float(2**(n-k))
    S_p = k*ghz_costs[mu_a] + (n-mu_a*k)
    R_par.append(R(S_p, D_100, kappa_opt(F_p)))


# ── Panel (b): n=100, k=1, vary mu ──────────────────────────────────────────

mu_vals = list(range(2, 21))
R_card_b, R_par_b = [], []
for mu in mu_vals:
    # Cardinality (k=1): |F| = mu * 2^(n-mu)  — decreases as mu grows
    F_c = float(mu) * float(2**(n-mu))
    S_c = float(dicke_costs[mu]) + (n-mu)
    R_card_b.append(R(S_c, D_100, kappa_opt(F_c)))

    # Parity (k=1): |F| = 2^(n-1)  — independent of mu
    F_p = float(2**(n-1))
    S_p = float(ghz_costs[mu]) + (n-mu)
    R_par_b.append(R(S_p, D_100, kappa_opt(F_p)))

# ── Standard (k=0) reference ────────────────────────────────────────────────
# uniform superposition, S=0, |F|=2^n
R_standard = R(100, D_100, kappa_opt(float(2**100)))


# ── Plot ─────────────────────────────────────────────────────────────────────

YLABEL = r'$S_\sigma + (D + 2S_\sigma)\,\kappa_\sigma^{\rm opt}$'
STD_LABEL = r'Standard'
STD_KW = dict(color='0.4', ls='-.', lw=1.8, zorder=1)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

# — Panel (a) —
ax1.axhline(R_standard, label=STD_LABEL, **STD_KW)
ax1.semilogy(k_vals, R_card, 'b-o', ms=6, lw=1.8, label='CA init. (cardinality)')
ax1.semilogy(k_vals, R_par,  'r-s', ms=6, lw=1.8, label='CA init. (parity)')
ax1.set_xlabel('Number of constraint sets $k$', fontsize=12)
ax1.set_ylabel(YLABEL, fontsize=11)
ax1.set_title(r'(a) $\mu=10$, $n=100$', fontsize=12)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_xticks(k_vals)

# — Panel (b) —
ax2.axhline(R_standard, label=STD_LABEL, **STD_KW)
ax2.semilogy(mu_vals, R_card_b, 'b-o', ms=6, lw=1.8, label='CA init. (cardinality)')
ax2.semilogy(mu_vals, R_par_b,  'r-s', ms=6, lw=1.8, label='CA init. (parity)')
ax2.set_xlabel('Constraint set size $\\mu$', fontsize=12)
ax2.set_ylabel(YLABEL, fontsize=11)
ax2.set_title(r'(b) $k=1$, $n=100$', fontsize=12)
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)
ax2.set_xticks(range(2, 21, 2))

plt.tight_layout()
plt.show()


