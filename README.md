# Grover's Algorithm with Constraint-Aware Initialization

Simulation code for the paper:

**Reducing Circuit Resources in Grover's Algorithm via Constraint-Aware Initialization**  
Eunok Bae, Jeonghyeon Shin, and Minjin Choi  
[arXiv:2601.17725](https://arxiv.org/abs/2601.17725)

---

## Requirements

```
pip install -r requirements.txt
```

Python 3.12, qiskit==2.1.1, qiskit-aer==0.17.1. 

---

## Structure

### `prob1/` — Cardinality constraints

Simulate Grover's algorithm on a 10-set exact-cover problem using cardinality-based constraint-aware initialization.

- `classical_preprocessing.py` — Preprocessing algorithm (Algorithm 1 in the paper)
- `prob1_std.py` — Standard (uniform) initialization
- `prob1_ca_c1.py` — CA initialization with one cardinality constraint (C1)
- `prob1_ca_c2.py` — CA initialization with a different single cardinality constraint (C2)
- `prob1_ca_c1c2.py` — CA initialization with two cardinality constraints (C1, C2)
- `prob1_ca_c1c2r3.py` — CA initialization with two cardinality constraints and one relaxed (overlap-bounded) constraint (C1, C2, R3)

### `prob2/` — Parity constraints (Figure 6)

Simulate Grover's algorithm on a weighted coverage problem using parity-based constraint-aware initialization.

- `classical_preprocessing.py` — Preprocessing algorithm (Algorithm 2 in the paper)
- `prob2_std.py` — Standard (uniform) initialization
- `prob2_c1.py` — CA initialization with one parity constraint (C1)
- `prob2_c1c2.py` — CA initialization with two parity constraints (C1, C2)

### `prob_extra/` — QFT-based oracle construction

Demonstrate that the oracle can be concretely implemented using a QFT-arithmetic approach (quantum adder), as an alternative to the ideal phase-flip oracle assumed in the paper. Applied to a smaller exact-cover instance (8 sets, 7 elements).

The oracle design follows:
- Quantum adder: A. Gilliam, S. Woerner, and C. Gonciulea, Quantum **5**, 428 (2021)
- Constraint-check: J.-R. Jiang and Y.-J. Wang, in *2023 VTS Asia Pacific Wireless Communications Symposium (APWCS)* (IEEE, 2023), pp. 1–5

- `classical_preprocessing.py` — Preprocessing for the smaller instance
- `prob_std.py` — Standard initialization with QFT oracle 
- `prob_c1c2.py` — CA initialization with QFT-based oracle 

### `resource_analysis.py` — Circuit resource analysis

Reproduce the resource analysis figure in the paper (gate count comparison across strategies as a function of constraint parameters). This figure does not appear in the current arXiv version and will be noted upon journal publication.
