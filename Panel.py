
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


def normalize_fluorophore(name):
    """Normalize fluorophore name to match similarity matrix."""
    name = name.strip()
    name = name.replace('/', '-')
    name = name.replace(' ', '')
    name = name.lower()
    return name


def fluor_to_similarity_name(name, sim_names_normalized, det_map=None):
    """Try to find a matching name in the similarity matrix.
    Returns (mapped_name, is_substituted) tuple."""
    normalized = normalize_fluorophore(name)
    if normalized in sim_names_normalized:
        return sim_names_normalized[normalized], False
    
    # Try detector-based substitution
    if det_map:
        for det, fl_list in det_map.items():
            for f in fl_list:
                if normalize_fluorophore(f) == normalized:
                    # Find a substitute from same detector that's in sim matrix
                    for candidate in fl_list:
                        candidate_norm = normalize_fluorophore(candidate)
                        if candidate_norm in sim_names_normalized:
                            return sim_names_normalized[candidate_norm], True
                    break
    
    return None, False


def load_similarity_matrix(sim_path):
    """Load similarity matrix and return DataFrame and name mapping."""
    import pandas as pd
    df = pd.read_csv(sim_path)
    df = df.set_index('Unnamed: 0')
    sim_names = df.index.tolist()
    sim_names_normalized = {normalize_fluorophore(n): n for n in sim_names}
    return df, sim_names


def load_detector_map(det_path):
    """Load detector-fluorophore mapping."""
    import pandas as pd
    df = pd.read_csv(det_path)
    det_map = {}
    for _, row in df.iterrows():
        det = row['Detector']
        fl = row['Fluorophore']
        if det not in det_map:
            det_map[det] = []
        det_map[det].append(fl)
    return det_map


def find_fluor_substitute(fluor, sim_names_normalized, det_map):
    """Find a substitute fluorophore that exists in the similarity matrix."""
    normalized = normalize_fluorophore(fluor)
    
    # Direct match
    if normalized in sim_names_normalized:
        return sim_names_normalized[normalized]
    
    # Try partial matching
    for norm_name, orig_name in sim_names_normalized.items():
        if normalized in norm_name or norm_name in normalized:
            return orig_name
    
    # Try detector-based substitution
    for det, fl_list in det_map.items():
        for f in fl_list:
            if normalize_fluorophore(f) == normalized:
                # Find a substitute from same detector that's in sim matrix
                for candidate in fl_list:
                    if candidate in sim_names_normalized:
                        return candidate
                break
    
    return fluor


def calculate_panel_complexity(fluorophores, sim_df, sim_names, det_map=None):
    """Calculate spectral complexity metrics for a panel."""
    import numpy as np
    
    sim_names_normalized = {normalize_fluorophore(n): n for n in sim_names}
    
    # Map fluorophores to similarity matrix names
    mapped_names = []
    unmapped_fluors = set()
    substitutions = {}
    for fl in fluorophores:
        mapped_fl, is_sub = fluor_to_similarity_name(fl, sim_names_normalized, det_map)
        if mapped_fl and mapped_fl in sim_df.columns:
            mapped_names.append(mapped_fl)
            if is_sub:
                substitutions[fl] = mapped_fl
        else:
            mapped_names.append(fl)
            unmapped_fluors.add(fl)
    
    n = len(mapped_names)
    
    if n < 2:
        return {
            'n_fluorophores': len(fluorophores),
            'matrix_complexity': 1.0,
            'unmapped_fluors': unmapped_fluors,
            'substitutions': substitutions
        }
    
    # Build full similarity matrix for all fluorophores
    S = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            fl1 = mapped_names[i]
            fl2 = mapped_names[j]
            
            if fl1 == fl2:
                S[i, j] = 1.0
            elif fl1 in sim_df.columns and fl2 in sim_df.columns:
                S[i, j] = sim_df.loc[fl1, fl2]
            elif fl1 in sim_df.columns and fl2 not in sim_df.columns:
                S[i, j] = 0.0
            elif fl1 not in sim_df.columns and fl2 in sim_df.columns:
                S[i, j] = 0.0
            else:
                S[i, j] = 0.0
    
    # Ensure symmetry
    S = (S + S.T) / 2
    np.fill_diagonal(S, 1.0)
    
    if not np.all(np.isfinite(S)):
        return {'matrix_complexity': np.inf, 'unmapped_fluors': unmapped_fluors, 'substitutions': substitutions}
    
    # A-optimality: mean variance inflation factor
    eigenvalues = np.linalg.eigvalsh(S)
    
    if not np.all(np.isfinite(eigenvalues)) or eigenvalues.min() <= 1e-12:
        matrix_complexity = np.inf
    else:
        matrix_complexity = np.mean(1.0 / eigenvalues)
    
    return {
        'n_fluorophores': len(fluorophores),
        'matrix_complexity': matrix_complexity,
        'unmapped_fluors': unmapped_fluors,
        'substitutions': substitutions
    }


def rank_combinations_by_complexity(assignment_result, sim_df, sim_names, det_map, top_n=10):
    """Rank combinations by spectral complexity (lowest first)."""
    sim_names_normalized = {normalize_fluorophore(n): n for n in sim_names}
    ranked = []
    
    for solution in assignment_result['combinations']:
        fluorophores = list(solution['assignments'].values())
        complexity = calculate_panel_complexity(fluorophores, sim_df, sim_names, det_map)
        
        ranked.append({
            'assignments': solution['assignments'],
            'unassigned': solution['unassigned'],
            'fluorophores': fluorophores,
            **complexity
        })
    
    ranked.sort(key=lambda x: x['matrix_complexity'])
    return ranked[:top_n], sim_names_normalized, det_map


# Example usage with CSV
import os
csv_path = os.path.join(os.path.dirname(__file__), "Data", "input.csv")
sim_path = os.path.join(os.path.dirname(__file__), "constants", "Similarity_Index.csv")
det_path = os.path.join(os.path.dirname(__file__), "constants", "Detector_Fluorophores.csv")

result = find_all_maximum_assignments(csv_path=csv_path, max_results=100000)

print(f"Maximum panel size: {result['max_assigned']}/{result['total_antibodies']}")
print(f"Total distinct combinations: {result['num_combinations']:,}")

# Load similarity data
sim_df, sim_names = load_similarity_matrix(sim_path)
det_map = load_detector_map(det_path)

# Rank by complexity
top_panels, sim_names_normalized, det_map = rank_combinations_by_complexity(result, sim_df, sim_names, det_map, top_n=10)

print("\n--- Top 10 Panels (Lowest Complexity) ---")
for i, panel in enumerate(top_panels, 1):
    print(f"\nRank {i} - A-optimality: {panel['matrix_complexity']:.4f}")
    for ab, fluor in panel['assignments'].items():
        tags = []
        if fluor in panel.get('substitutions', {}):
            sub = panel['substitutions'][fluor]
            tags.append(f"[substituted: {sub}]")
        if fluor in panel.get('unmapped_fluors', set()):
            tags.append("[unmapped]")
        tag = " " + " ".join(tags) if tags else ""
        print(f"  {ab:12} -> {fluor}{tag}")
