from collections import defaultdict

def build_index_sets(list_of_lists):
    v2idx = defaultdict(set)
    for i, lst in enumerate(list_of_lists):
        for v in lst:
            v2idx[v].add(i+1)
    values = sorted(v2idx.keys())
    constraints = [v2idx[v] for v in values]
    return values, constraints

# if __name__ == "__main__":
#     A = [
#         [1, 5],  # A1
#         [1, 3, 6],  # A2
#         [1, 2, 5],  # A3
#         [1, 7],  # A4
#         [3, 4],  # A5
#         [4, 6],  # A6
#         [2, 4, 5],  # A7
#         [2, 7],  # A8
#         [6, 7],  # A9
#         [3, 5, 7],  # A10
#     ]
#     values, constraints = build_index_sets(A)
#     for idx, subset in zip(values, constraints):
#         print(f"value {idx}: {subset}")

from math import comb
def greedy_search(constraints, cardinalities, tau=2):
    ordered = []
    for i, C in enumerate(constraints):
        Ci = set(C)
        k = len(Ci)
        b = int(cardinalities[i])
        r = float('inf') if k == 0 else comb(k, b) / 2**k
        ordered.append((i, Ci, b, r))
    ordered.sort(key=lambda t: t[3])

    used = set()
    organized = []

    # Phase 1: disjoint subset
    for idx, subset, b, _ in ordered:
        if subset.isdisjoint(used):
            organized.append((idx, set(subset), f"b == {b}"))
            used.update(subset)

    # Phase 2: bounded-overlap
    chosen_ids = {idx for idx, _, _ in organized}
    remaining = [(i, Ci, b) for (i, Ci, b, _) in ordered if i not in chosen_ids]
    for idx, subset, b in remaining:
        overlap = subset & used
        if len(overlap) <= tau:
            adjusted = subset - used
            if adjusted:
                lb = max(0, b-len(overlap))
                ub = b
                organized.append((idx, adjusted, f"{lb} ≤ b ≤ {ub}"))
                used.update(adjusted)

    return organized

if __name__ == "__main__":
    A = [
        [1, 5],  # A1
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
    values, constraints = build_index_sets(A)
    b = [1 for _ in range(len(values))]
    for idx, subset in zip(values, constraints):
        print(f"value {idx}: {subset}")

    organized = greedy_search(constraints=constraints, cardinalities=b, tau=1)

    print("=== Organized sets ===")
    for idx, subset, cond in organized:
        print(f"C{idx+1}: {subset},  condition: {cond}")

