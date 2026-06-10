import csv

max_energy = float('-inf')
min_energy = float('inf')
max_hop = None
min_hop = None

with open('final_results.csv', mode='r') as file:
    # Read the CSV file as a dictionary
    reader = csv.DictReader(file)
    
    for row in reader:
        hop = int(row['hop'])
        energy = float(row['energy_eV'])
        
        # Check for maximum energy
        if energy > max_energy:
            max_energy = energy
            max_hop = hop
            
        # Check for minimum energy
        if energy < min_energy:
            min_energy = energy
            min_hop = hop

print(f"Hop {max_hop}: {max_energy} (Maximum Energy)")
print(f"Hop {min_hop}: {min_energy} (Minimum Energy)")
