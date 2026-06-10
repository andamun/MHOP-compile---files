import matplotlib.pyplot as plt
import numpy as np
import io
import os
import re

# import data
LOG_FILE_NAME = 'md.log'
T_target = int(input('Temperature_target:'))

try:
     with open(LOG_FILE_NAME, 'r') as f:
         log_data = f.read()
except FileNotFoundError:
     print(f"Error: Log file '{LOG_FILE_NAME}' not found.")
     exit()

# --- Data Loading and Cleaning ---
# Use io.StringIO to treat the string data as a file object
# The data is cleaned by replacing multiple spaces/tabs with single spaces
# to ensure loadtxt reads it correctly.
data_clean = re.sub(r'[ \t]+', ' ', log_data.strip())
data_file = io.StringIO(data_clean)

# Load data using numpy's loadtxt.
# skiprows=1 skips the column headers.
# unpack=True loads columns into separate variables.
try:
    Time, Etot, Epot, Ekin, Temp = np.loadtxt(
        data_file,
        skiprows=1,
        unpack=True
    )
except ValueError as e:
    # If using the file in a real environment, this error handling is critical
    print(f"Error loading data from {LOG_FILE_NAME}. Check for non-numeric values or inconsistent columns: {e}")
    exit()

# --- Plotting Function Definitions ---

# 1. Plot Time vs. Etot
def plot_etot_vs_time(Time, energy, energy_name):
    """Generates and displays the Total Energy vs. Time plot."""
    plt.figure(figsize=(10, 6))
    plt.plot(Time, energy, label=f'{energy_name}', color='darkred', linewidth=2)
    plt.xlabel('Time (ps)', fontsize=14)
    plt.ylabel('Total Energy (eV)', fontsize=14)
    plt.title(f'Molecular Dynamics: Total Energy {energy_name} vs. Time', fontsize=16)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig(f'plot_{energy_name}.png')
    plt.show()


# 2. Plot Time vs. All 3 Energies
def plot_all_energies_vs_time(Time, Etot, Epot, Ekin):
    """Generates and displays the plot of all three energy components vs. Time."""
    plt.figure(figsize=(10, 6))

    # Plot Epot, typically the most important for stability
    plt.plot(Time, Epot, label='Potential Energy ($E_{pot}$)', color='blue', linewidth=1.5, alpha=0.8)

    # Plot Ekin, showing fluctuations due to temperature
    plt.plot(Time, Ekin, label='Kinetic Energy ($E_{kin}$)', color='green', linewidth=1.5, alpha=0.8)

    # Plot Etot, checking for drift (should be nearly constant in NVT)
    plt.plot(Time, Etot, label='Total Energy ($E_{tot}$)', color='darkred', linewidth=2, linestyle='--')

    plt.xlabel('Time (ps)', fontsize=14)
    plt.ylabel('Energy (eV)', fontsize=14)
    plt.title('Molecular Dynamics: Energy Components vs. Time', fontsize=16)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(fontsize=12, loc='best')
    plt.tight_layout()
    plt.savefig('plot_2_all_energies.png')
    plt.show()

# 3. Plot Time vs. Temperature
def plot_temp_vs_time(Time, Temp,T_target):
    """Generates and displays the Temperature vs. Time plot."""
    plt.figure(figsize=(10, 6))
    T_target = T_target

    # Plot the calculated temperature
    plt.plot(Time, Temp, label='Instantaneous Temperature', color='orange', linewidth=2)

    # Plot a horizontal line for the target temperature (NVT ensemble)
    plt.axhline(T_target, color='gray', linestyle='--', label=f'Target Temperature ({T_target} K)')

    plt.xlabel('Time (ps)', fontsize=14)
    plt.ylabel('Temperature (K)', fontsize=14)
    plt.title('Molecular Dynamics: Temperature Control (NVT)', fontsize=16)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig('plot_3_temperature.png')
    plt.show()


# --- Execution ---
if __name__ == '__main__':
    print("Generating Plot 1: Total Energy vs. Time...")
    plot_etot_vs_time(Time, Etot,'Total Energy')
    plot_etot_vs_time(Time, Epot,'Potential Energy')
    plot_etot_vs_time(Time, Ekin,'Kinetic Energy')
    
    print("Generating Plot 2: All Energy Components vs. Time...")
    plot_all_energies_vs_time(Time, Etot, Epot, Ekin)
    
    print("Generating Plot 3: Temperature vs. Time...")
    plot_temp_vs_time(Time, Temp,T_target)

    print("\nPlots saved as 'plot_1_etot.png', 'plot_2_all_energies.png', and 'plot_3_temperature.png'.")


