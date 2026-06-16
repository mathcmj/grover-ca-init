from collections import defaultdict

def build_index_sets(list_of_lists, even_coeff_list=None):
    v2idx = defaultdict(set)
    for i, lst in enumerate(list_of_lists):
        if even_coeff_list is not None and i in even_coeff_list:
            continue
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
#     values, constraints = build_index_sets(A, [2, 5, 9])
#     for idx, subset in zip(values, constraints):
#         print(f"value {idx}: {subset}")

def greedy_search(constraints):
    ordered = []
    for i, C in enumerate(constraints):
        Ci = set(C)
        k = len(Ci)
        ordered.append((i, Ci, k))
    ordered.sort(key=lambda t: t[2])

    used = set()
    organized = []

    # Phase 1: disjoint subset
    for idx, subset, _ in ordered:
        if subset.isdisjoint(used):
            organized.append((idx, set(subset)))
            used.update(subset)
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
    values, constraints = build_index_sets(A, [2, 5, 9])
    b = [1 for _ in range(len(values))]
    for idx, subset in zip(values, constraints):
        print(f"value {idx}: {subset}")

    organized = greedy_search(constraints=constraints)

    print("=== Organized sets ===")
    for idx, subset in organized:
        print(f"C{idx+1}: {subset}")

