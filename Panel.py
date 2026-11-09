
def assign_fluorophores(antibodies=None, csv_path=None):
    """
    Assign fluorophores to antibodies using greedy matching.
    
    Args:
        antibodies: dict mapping antibody names to list of fluorophores
        csv_path: path to CSV file with 'Antibody' and 'Fluorophore' columns
    
    If csv_path is provided, loads antibody-fluorophore pairs from the CSV.
    """
    import os
    if csv_path and not antibodies:
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"CSV not found: {csv_path}")
        antibodies = {}
        import csv
        with open(csv_path, 'r', newline='') as f:
            reader = csv.DictReader(f)
            for row in reader:
                ab = row['Antibody'].strip()
                fl = row['Fluorophore'].strip()
                if ab not in antibodies:
                    antibodies[ab] = []
                if fl not in antibodies[ab]:
                    antibodies[ab].append(fl)

    fluor_to_antibody = {}

    def find_assignment(antibody, visited):
        for fluor in antibodies[antibody]:

            if fluor in visited:
                continue

            visited.add(fluor)

            if (fluor not in fluor_to_antibody or
                find_assignment(fluor_to_antibody[fluor], visited)):

                fluor_to_antibody[fluor] = antibody
                return True

        return False

    # Attempt to assign every antibody
    for antibody in antibodies:
        find_assignment(antibody, set())

    # Reverse mapping
    assignments = {
        antibody: fluor
        for fluor, antibody in fluor_to_antibody.items()
    }

    unassigned = [
        antibody for antibody in antibodies
        if antibody not in assignments
    ]

    return assignments, unassigned


def find_all_maximum_assignments(antibodies=None, csv_path=None, max_results=100):
    """
    Find all possible assignment combinations that achieve the maximum number of assignments.
    
    Args:
        antibodies: dict mapping antibody names to list of fluorophores
        csv_path: path to CSV file with 'Antibody' and 'Fluorophore' columns
        max_results: maximum number of combinations to return
    
    Returns:
        dict with:
            - max_assigned: maximum number of antibodies that can be assigned
            - total_antibodies: total number of antibodies
            - num_combinations: total number of distinct maximum assignments
            - combinations: list of solution dicts (up to max_results)
            - truncated: whether the solution list was truncated
    """
    import os
    
    if csv_path and not antibodies:
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"CSV not found: {csv_path}")
        antibodies = {}
        import csv
        with open(csv_path, 'r', newline='') as f:
            reader = csv.DictReader(f)
            for row in reader:
                ab = row['Antibody'].strip()
                fl = row['Fluorophore'].strip()
                if ab not in antibodies:
                    antibodies[ab] = []
                if fl not in antibodies[ab]:
                    antibodies[ab].append(fl)

    # Deduplicate fluorophores while preserving order
    options = {ab: list(dict.fromkeys(fls)) for ab, fls in antibodies.items()}
    original_order = list(antibodies)

    # Sort by most-constrained-first for efficient search
    search_order = sorted(
        original_order,
        key=lambda ab: (len(options[ab]), original_order.index(ab))
    )

    # Find maximum matching size using greedy
    fluor_to_antibody = {}

    def find_assignment(antibody, visited):
        for fluor in options[antibody]:
            if fluor in visited:
                continue
            visited.add(fluor)
            if fluor not in fluor_to_antibody or find_assignment(fluor_to_antibody[fluor], visited):
                fluor_to_antibody[fluor] = antibody
                return True
        return False

    for ab in search_order:
        find_assignment(ab, set())

    max_count = len(fluor_to_antibody)

    # Collect all fluorophores
    all_fluors = set()
    for fl_list in options.values():
        all_fluors.update(fl_list)
    all_fluors = sorted(all_fluors)
    fluor_idx = {fl: i for i, fl in enumerate(all_fluors)}

    # Fast counting using bitmask DP (original order, no sorting)
    from functools import lru_cache

    ab_list_for_count = list(options.keys())
    count_total = len(ab_list_for_count)
    count_masks = []
    for ab in ab_list_for_count:
        mask = 0
        for fl in options[ab]:
            mask |= (1 << fluor_idx[fl])
        count_masks.append(mask)

    @lru_cache(maxsize=None)
    def count_ways(idx, used_mask, target):
        if target == 0:
            return 1
        remaining = count_total - idx
        if remaining < target:
            return 0
        if idx == count_total:
            return 0
        total_count = 0
        mask = count_masks[idx]
        available = mask & ~used_mask
        while available:
            fl_bit = available & (-available)
            total_count += count_ways(idx + 1, used_mask | fl_bit, target - 1)
            available &= available - 1
        total_count += count_ways(idx + 1, used_mask, target)
        return total_count

    total_combinations = count_ways(0, 0, max_count)

    # Enumerate all maximum assignments via backtracking (sorted order)
    combinations = []
    truncated = False

    def enumerate_ways(idx, assigned, used_fluors):
        nonlocal truncated
        remaining = len(search_order) - idx

        if len(assigned) + remaining < max_count:
            return

        if len(assigned) == max_count:
            assignments = {ab: assigned[ab] for ab in original_order if ab in assigned}
            unused = sorted(fl for fl in all_fluors if fl not in used_fluors)
            unassigned = [ab for ab in original_order if ab not in assigned]
            combinations.append({
                'assignments': assignments,
                'unassigned': unassigned,
                'unused_fluorophores': unused
            })
            if len(combinations) >= max_results:
                truncated = True
            return

        if idx == len(search_order):
            return

        ab = search_order[idx]

        for fluor in options[ab]:
            if fluor in used_fluors:
                continue
            assigned[ab] = fluor
            used_fluors.add(fluor)
            enumerate_ways(idx + 1, assigned, used_fluors)
            del assigned[ab]
            used_fluors.remove(fluor)
            if truncated:
                return

        enumerate_ways(idx + 1, assigned, used_fluors)

    enumerate_ways(0, {}, set())

    return {
        'max_assigned': max_count,
        'total_antibodies': len(antibodies),
        'num_combinations': total_combinations,
        'combinations': combinations[:max_results],
        'truncated': truncated
    }


# Example usage with CSV
import os
csv_path = os.path.join(os.path.dirname(__file__), "Data", "antibody_fluorophores.csv")

result = find_all_maximum_assignments(csv_path=csv_path)

print(f"Maximum panel size: {result['max_assigned']}/{result['total_antibodies']}")
print(f"Total distinct combinations: {result['num_combinations']:,}")

print("\n--- Example Combination ---")
for ab, fluor in result['combinations'][0]['assignments'].items():
    print(f"{ab:12} -> {fluor}")
