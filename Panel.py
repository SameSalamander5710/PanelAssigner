
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
            - max_count: maximum number of antibodies that can be assigned
            - total_antibodies: total number of antibodies
            - num_combinations: total number of distinct maximum assignments
            - combinations: list of assignment dicts (up to max_results)
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

    ab_list = list(antibodies.keys())
    total = len(ab_list)

    # Find maximum matching size using greedy
    def greedy_max():
        fluor_to_antibody = {}
        def find_assignment(antibody, visited):
            for fluor in antibodies[antibody]:
                if fluor in visited:
                    continue
                visited.add(fluor)
                if fluor not in fluor_to_antibody or find_assignment(fluor_to_antibody[fluor], visited):
                    fluor_to_antibody[fluor] = antibody
                    return True
            return False
        for antibody in ab_list:
            find_assignment(antibody, set())
        return len(fluor_to_antibody)

    max_count = greedy_max()

    # Collect all fluorophores
    all_fluors = set()
    for fl_list in antibodies.values():
        all_fluors.update(fl_list)
    all_fluors = sorted(all_fluors)
    fluor_idx = {fl: i for i, fl in enumerate(all_fluors)}

    # Precompute antibody fluor masks
    ab_masks = []
    for ab in ab_list:
        mask = 0
        for fl in antibodies[ab]:
            mask |= (1 << fluor_idx[fl])
        ab_masks.append(mask)

    # Memoized counting: count ways to assign k antibodies from idx..total-1
    # using only fluorophores not in used_mask
    from functools import lru_cache

    @lru_cache(maxsize=None)
    def count_ways(idx, used_mask, target):
        """Count ways to assign exactly 'target' antibodies from idx onwards."""
        if target == 0:
            return 1
        remaining = total - idx
        if remaining < target:
            return 0
        if idx == total:
            return 0
        
        total_count = 0
        # Try assigning antibody[idx] to each available fluor
        mask = ab_masks[idx]
        available = mask & ~used_mask
        while available:
            fl_bit = available & (-available)  # lowest set bit
            total_count += count_ways(idx + 1, used_mask | fl_bit, target - 1)
            available &= available - 1  # remove lowest set bit
        
        # Try skipping antibody[idx]
        total_count += count_ways(idx + 1, used_mask, target)
        
        return total_count

    total_combinations = count_ways(0, 0, max_count)

    # Enumerate combinations
    combinations = []

    def enumerate_ways(idx, assigned, used_mask, remaining):
        nonlocal total_combinations
        if len(combinations) >= max_results:
            return
        
        if remaining == 0:
            combinations.append(dict(assigned))
            return
        
        if idx == total:
            return
        
        # Prune: check if enough ways remain
        ways = count_ways(idx, used_mask, remaining)
        if ways == 0:
            return

        # Try assigning antibody[idx] to each available fluor
        mask = ab_masks[idx]
        available = mask & ~used_mask
        while available:
            fl_bit = available & (-available)
            fl = all_fluors[fl_bit.bit_length() - 1]
            assigned[ab_list[idx]] = fl
            enumerate_ways(idx + 1, assigned, used_mask | fl_bit, remaining - 1)
            del assigned[ab_list[idx]]
            if len(combinations) >= max_results:
                return
            available &= available - 1
        
        # Try skipping antibody[idx]
        enumerate_ways(idx + 1, assigned, used_mask, remaining)

    enumerate_ways(0, {}, 0, max_count)

    return {
        'max_count': max_count,
        'total_antibodies': total,
        'num_combinations': total_combinations,
        'combinations': combinations[:max_results]
    }


# Example usage with CSV
import os
csv_path = os.path.join(os.path.dirname(__file__), "Data", "antibody_fluorophores.csv")

result = find_all_maximum_assignments(csv_path=csv_path)

print(f"Maximum assignments: {result['max_count']}/{result['total_antibodies']}")
num = result['num_combinations']
if num > 999999:
    print(f"Total distinct combinations: >999,999")
else:
    print(f"Total distinct combinations: {num:,}")
print(f"Combinations shown: {len(result['combinations'])}")

print("\n--- Combination 1 ---")
for ab, fl in result['combinations'][0].items():
    print(f"{ab:10} -> {fl}")
