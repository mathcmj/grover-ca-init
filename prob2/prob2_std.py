import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_aer import Aer

n_s = 10

qc_Grover = QuantumCircuit(n_s, n_s)

set_collection = [[1, 5],  # A1
                  [1, 3, 6],  # A2
                  [1, 2, 5],  # A3
                  [1, 7],  # A4
                  [3, 4],  # A5
                  [4, 6],  # A6
                  [2, 4, 5],  # A7
                  [2, 7],  # A8
                  [6, 7],  # A9
                  [3, 5, 7],  # A10
                 ]


from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import Decompose

def mcx_noancilla_as_instrs(num_ctrls: int, basis=('rz', 'h','cx')):

    if num_ctrls < 1:
        raise ValueError("num_ctrls must be >= 1")
    k = num_ctrls

    sub = QuantumCircuit(k + 1)
    sub.mcx(list(range(k)), k)

    sub_u = transpile(sub, basis_gates=list(basis), optimization_level=0)

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

from qiskit.circuit.library import UnitaryGate
def oracle():
    solutions = ['0111010011']
    diag = np.ones(2**n_s, dtype=complex)
    for sol in solutions:
        if isinstance(sol, str):
            idx = int(sol, 2)
        else:
            idx = int(''.join(str(b) for b in sol), 2)
        diag[idx] = -1.0
    U = np.diag(diag)
    qc_Grover.append(UnitaryGate(U, label='Oracle'), [qc_Grover.qubits[i] for i in range(n_s)])


def reflection():
    for i in range(n_s-1):
        qc_Grover.h(i)
        qc_Grover.x(i)
    qc_Grover.z(n_s-1)
    transpiled_mcx(qc_Grover, [i for i in range(n_s-1)], n_s-1)
    # qc_Grover.mcx(control_qubits=[i for i in range(n_s-1)], target_qubit=n_s-1)
    qc_Grover.z(n_s-1)
    for i in range(1, n_s-1):
        qc_Grover.x(i)
    qc_Grover.z(0)
    qc_Grover.x(0)
    qc_Grover.z(0)
    for i in range(n_s-1):
        qc_Grover.h(i)
    qc_Grover.barrier()

## Grover algorithm
for i in range(n_s):
    qc_Grover.h(i)
qc_Grover.barrier()
# iter = 25
iter = round((np.pi/(4*np.arcsin(np.sqrt(1/2**n_s))))-1/2)
print(iter)
for i in range(iter):
    oracle()
    reflection()

qc_Grover.measure(list(range(n_s)), list(reversed(range(n_s))))

# ## TEST
# # print(qc_Grover.draw(output='text'))
# n_shots = 1000
# backend = Aer.get_backend('aer_simulator')
# job = backend.run(qc_Grover, shots=n_shots)
# result = job.result()
# counts = result.get_counts()
# sorted_counts = dict(sorted(counts.items(), key=lambda item: item[1], reverse=True))
# sol_outcome = '1100101110'
# sol_count = sorted_counts.get(sol_outcome, 0)
# print('results:', sorted_counts)
# print('sol_count:', sol_count)
# print("gate_type: ", qc_Grover.count_ops())


## emulator (with depolarizing noise)
from qiskit_aer.noise import NoiseModel, depolarizing_error

sq_error = 0.00001
tq_error = 0.0001
backend = Aer.get_backend('aer_simulator')
noise_model = NoiseModel()
single_qubit_error = depolarizing_error(sq_error, 1)
two_qubit_error = depolarizing_error(tq_error, 2)
noise_model.add_all_qubit_quantum_error(single_qubit_error, ['rz', 'z', 'h', 'x'])
noise_model.add_all_qubit_quantum_error(two_qubit_error, ['cx', 'cry'])


f_ideal = open('results_prob2_std_ideal', 'a')
f_noisy = open('results_prob2_std_noisy', 'a')

print('iter:', iter, file=f_ideal)
print('iter:', iter, file=f_noisy)
print('sq_error:', sq_error, 'tq_error:', tq_error, file=f_noisy)
sol_count_ideal = []
sol_count_noisy = []
for _ in range(5):
    n_shots = 1000
    job_ideal = backend.run(qc_Grover, shots=n_shots)
    job_noise = backend.run(qc_Grover, noise_model=noise_model, shots=n_shots)
    result_ideal = job_ideal.result()
    result_noise = job_noise.result()
    counts_ideal = result_ideal.get_counts()
    counts_noise = result_noise.get_counts()
    sorted_counts_ideal = dict(sorted(counts_ideal.items(), key=lambda item: item[1], reverse=True))
    sorted_counts_noise = dict(sorted(counts_noise.items(), key=lambda item: item[1], reverse=True))
    print('results(ideal):', sorted_counts_ideal, file=f_ideal)
    print('results(noise):', sorted_counts_noise, file=f_noisy)
    sol_outcome = '1100101110'
    sol_count_ideal.append(sorted_counts_ideal.get(sol_outcome, 0))
    sol_count_noisy.append(sorted_counts_noise.get(sol_outcome, 0))
print('sol_count_ideal:', sol_count_ideal, file=f_ideal)
print('sol_count_noisy:', sol_count_noisy, file=f_noisy)
print('='*100, file=f_ideal)
print('='*100, file=f_noisy)