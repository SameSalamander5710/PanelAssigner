
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


# Example usage with CSV
import os
csv_path = os.path.join(os.path.dirname(__file__), "Data", "antibody_fluorophores.csv")
assignments, unassigned = assign_fluorophores(csv_path=csv_path)
total = len(assignments) + len(unassigned)

print("ASSIGNMENTS")
for antibody, fluor in assignments.items():
    print(f"{antibody:10} -> {fluor}")

print("\nUNASSIGNED")
print(unassigned)

print(f"\nAssigned: {len(assignments)}/{total}")
