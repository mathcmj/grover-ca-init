import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_aer import Aer

n_s = 8
n_e = 5
n_u = 2
n_q = n_s + n_e*n_u + 1    # n_q = 19

qc_Grover = QuantumCircuit(n_q, n_s)

set_collection = [[1, 5],  # A1
                  [1, 3, 6],  # A2
                  [1, 2, 7],  # A3
                  [2, 4],  # A4
                  [3, 4],  # A5
                  [3, 4, 6],  # A6
                  [6, 7],  # A7
                  [5, 7],  # A8
                  ]

U = [2, 3, 5, 6, 7]
for subset in set_collection:
    if 1 in subset:
        subset.remove(1)
    if 4 in subset:
        subset.remove(4)

mapping = {v: i+1 for i, v in enumerate(U)}
set_collection = [[mapping[x] for x in subset] for subset in set_collection]
# print(set_collection)


def oracle():
    for idx in range(len(set_collection)):
        lst = set_collection[idx]
        for i in lst:
            for j in range(n_u):
                qc_Grover.cp(theta=2*np.pi*(1/2**(j+1)), control_qubit=[idx], target_qubit=n_s+n_u*(i-1)+j)
            qc_Grover.barrier()
    ## QFT^{dagger}
    for i in range(n_e):
        for j in range(n_u):
            for k in range(j):
                qc_Grover.cp(theta=-2*np.pi*(1/2**(j+1-k)), control_qubit=n_s+n_u*i+k, target_qubit=n_s+n_u*i+j)
            qc_Grover.h(n_s+n_u*i+j)
    qc_Grover.barrier()
    ## minus sign
    for i in range(n_e):
        for j in range(n_u-1):
            qc_Grover.x(n_s+n_u*(i+1)-1-j)
    qc_Grover.mcx(control_qubits=[n_s+i for i in range(n_e*n_u)], target_qubit=n_q-1)
    for i in range(n_e):
        for j in range(n_u-1):
            qc_Grover.x(n_s+n_u*(i+1)-1-j)
    qc_Grover.barrier()
    ## QFT
    for i in range(n_e):
        for j in reversed(range(n_u)):
            qc_Grover.h(n_s+n_u*i+j)
            for k in reversed(range(j)):
                qc_Grover.cp(theta=2*np.pi*(1/2**(j+1-k)), control_qubit=n_s+n_u*i+k, target_qubit=n_s+n_u*i+j)
    qc_Grover.barrier()
    for idx in range(len(set_collection)):
        lst = set_collection[idx]
        for i in lst:
            for j in reversed(range(n_u)):
                qc_Grover.cp(theta=-2*np.pi*(1/2**(j+1)), control_qubit=[idx], target_qubit=n_s+n_u*(i-1)+j)
            qc_Grover.barrier()

def scs_n_1(circ, m, i):
    circ.cry(theta=2*np.arccos(np.sqrt(1/m)), control_qubit=i+1, target_qubit=i)
    circ.cx(control_qubit=i, target_qubit=i+1)
def inv_scs_n_1(circ, m, i):
    circ.cx(control_qubit=i, target_qubit=i+1)
    circ.cry(theta=-2*np.arccos(np.sqrt(1/m)), control_qubit=i+1, target_qubit=i)

def unitary_W(circ, q_list):
    n = len(q_list)
    circ.x(q_list[-1])
    for iter in range(n-1):
        scs_n_1(circ=circ, m=n-iter, i=q_list[n-2-iter])
def inv_unitary_W(circ, q_list):
    n = len(q_list)
    for iter in reversed(range(n-1)):
        inv_scs_n_1(circ=circ, m=n-iter, i=q_list[n-2-iter])
    circ.x(q_list[-1])


def reflection():
    inv_unitary_W(qc_Grover, [0, 1, 2])
    inv_unitary_W(qc_Grover, [3, 4, 5])
    qc_Grover.h(7)
    qc_Grover.h(6)
    qc_Grover.barrier()
    for i in range(n_s):
        qc_Grover.x(i)
    qc_Grover.h(n_s-1)
    qc_Grover.mcx(control_qubits=[i for i in range(n_s-1)], target_qubit=n_s-1)
    qc_Grover.h(n_s-1)
    for i in range(1, n_s):
        qc_Grover.x(i)
    qc_Grover.z(0)
    qc_Grover.x(0)
    qc_Grover.z(0)
    qc_Grover.barrier()
    unitary_W(qc_Grover, [0, 1, 2])
    unitary_W(qc_Grover, [3, 4, 5])
    qc_Grover.h(7)
    qc_Grover.h(6)
    qc_Grover.barrier()

## Grover algorithm
unitary_W(qc_Grover, [0, 1, 2])
unitary_W(qc_Grover, [3, 4, 5])
qc_Grover.h(7)
qc_Grover.h(6)
for i in range(n_e*n_u):
    qc_Grover.h(n_s+i)
qc_Grover.x(n_q-1)
qc_Grover.h(n_q-1)
qc_Grover.barrier()
# iter = 4
iter = round((np.pi/(4*np.arcsin(np.sqrt(1/(3*3*2*2)))))-1/2)
# print(iter)
for i in range(iter):
    oracle()
    reflection()

qc_Grover.measure(list(range(n_s)), list(reversed(range(n_s))))

## TEST
# print(qc_Grover.draw(output='text'))
n_shots = 1000
backend = Aer.get_backend('aer_simulator')
job = backend.run(qc_Grover, shots=n_shots)
result = job.result()
counts = result.get_counts()
sorted_counts = dict(sorted(counts.items(), key=lambda item: item[1], reverse=True))
sol_outcome = '01010001'
sol_count = sorted_counts.get(sol_outcome, 0)
print('results:', sorted_counts)
print('sol_count:', sol_count)
print("gate_type: ", qc_Grover.count_ops())

