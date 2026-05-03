### CALPHAD Practice
# Mitchell Shapiro-Albert
# May 1, 2026

# _________ Imports _________
print("Beginning code...")

# Handling file imports
from pathlib import Path
# For graphing
import matplotlib.pyplot as plt
# Main package to be used
from pycalphad import Database, binplot, ternplot, equilibrium, calculate
import pycalphad.variables as v # Import specific submodules as needed
from pycalphad.plot.utils import phase_legend
# To remove distracting warning messages; see note below
import warnings
# To allow for timing functions/processes
import time
import threading
# For math
import numpy as np
import xarray as xr

print("Imports successful.")

# _________ Warning Handling _________
# NOTE: you should remove the warning handling if you are using this for research purposes! These suppress warning
# messages that you should be aware of as a user. E.g., if you are calculating molar volume or magnetic heat capacity,
# these warnings will tell you that the .tdb file you are using is not a good candidate.
warnings.filterwarnings("ignore", message=".*Type definitions using IF/THEN logic.*")
warnings.filterwarnings("ignore", message=".*no corresponding TYPE_DEFINITION line was found.*")

# _________ Functions _________
print("Writing functions...")

# Write a function that can run and periodically print a message to show users the code is not frozen
def status_update(stop_event, interval=5):
    start_time = time.time()
    while not stop_event.is_set():
        time.sleep(interval)
        if not stop_event.is_set():
            elapsed = time.time() - start_time
            print(f"... still cooking over here: ({elapsed:.0f}s elapsed) ...")

# This is called by adding the following "wrapper" to any function:

# start_time = time.time()
# print("Start message")
# stop_heartbeat = threading.Event()
# heartbeat_thread = threading.Thread(target=print_status, args=(stop_heartbeat, 5)) # Prints every 5s
# heartbeat_thread.start()
# \\\ Entire function in a "try:" block \\\
# finally:
#   stop_heartbeat.set()
#   heartbeat_thread.join()
#   end_time = time.time()
#   print("End message.")

# Write a function to retrieve the .tbd file from the subdirectory
def retrieve_tdb(file_name):
    script_dir = Path(__file__).resolve()
    target_path = script_dir.parent / "TDB_Files" / f"{file_name}"
    if target_path.exists():
        return target_path
    else:
        raise FileNotFoundError(f"{target_path} does not exist")

# Write a function to fully summarize a Database object
def orient_database(dbf, sum_length="Full"):
    # Print header
    print(f"{'=' * 20} IMPORTED DATABASE SUMMARY {'=' * 20}")
    # Print basic components and species
    print(f"\n[COMPONENTS & SPECIES]")
    print(f"Elements: {sorted(dbf.elements)}")
    print(f"Species:  {[s.name for s in dbf.species]}")
    # If the user wants a full preview of the database loaded:
    if sum_length == "Full":
        # Print parameters and functions (Symbols); includes Gibbs energy, Curie temperature, mathematical relations
        print(f"\n[SYMBOLS/FUNCTIONS]")
        functions = [sym for sym in dbf.symbols.keys()]
        print(f"Total Symbols: {len(functions)}")
        if functions:
            print(f"Sample Symbols: {functions[:5]}...")
    # Print information for each relevant phase
    print(f"\n[PHASE DESCRIPTIONS]")
    phase_names = sorted(dbf.phases.keys())
    print(f"Total Phases: {len(phase_names)}")
    # Loop through each relevant phase:
    for name in phase_names:
        phase = dbf.phases[name]
        # Extract sublattice and constituent atom information
        sublattices = phase.sublattices
        # Constituents are stored as a tuple of species for each sublattice
        constituents = [[s.name for s in c] for c in phase.constituents]
        # Print the name of the phase to terminal
        print(f"\nPhase: {name}")
        # If calling full summary, print the above findings for each phase to terminal
        if sum_length == "Full":
            print(f"  - Sublattices: {sublattices}")
            for i, (sites, species) in enumerate(zip(sublattices, constituents)):
                print(f"    Sublattice {i + 1} ({sites} sites): {', '.join(species)}")
            # Check for model hints (e.g., if it's a magnetic or ordered phase)
            if hasattr(phase, 'model_hints') and phase.model_hints:
                print(f"    Hints: {list(phase.model_hints.keys())}")
    # Neatly format the end of the output
    print(f"\n{'=' * 60}")

# Write a function to plot a binary phase diagram given a retrieved_tdb object, components, and physical conditions:
# Ex: components = ["AL", "NI", "VA"], x_component = "AL"
def plot_binary_diagram(dbf, components, x_component, x_step_range=(0, 1, 0.05), phases=None, moles=1, p_pa=101325, temp_range=(300, 2000, 10), save=False):
    print("Graphing...")
    # If phases are not provided, gets all possible phases programmatically with minor type error handling
    if phases is None or not isinstance(phases, list):
        phases = list(dbf.phases.keys())
    # Create axes object using pycalphad built-in binplot
    ax = binplot(dbf, components, phases, {v.N: moles, v.P: p_pa, v.T: temp_range, v.X(x_component): x_step_range},
            plot_kwargs={"tielines": True}) # Keeps ugly "shading" in graph region
    # Add some flashy formatting for the graph
    ax.set_title(f"Binary Phase Diagram: {'-'.join([c for c in components if c != 'VA'])}", fontsize=14, pad=15)
    ax.set_xlabel(f"Mole Fraction of {x_component}", fontsize=12)
    ax.set_ylabel("Temperature (K)", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.6) # Add a subtle grid for readability
    if save:
        plt.savefig("binary_diagram.jpg", format="jpg", dpi=300, bbox_inches="tight")
        print("Binary phase diagram saved as binary_diagram.jpg.")
    plt.show()

# Write a function to plot a ternary phase diagram given a retrieved_tdb object, components, and physical conditions:
def plot_ternary_diagram(dbf, step_range=(0, 1, 0.05), temp=1000, phases=None, p_pa=101325, save=False):
    print("Graphing...")
    # Get the pure elements for labels
    pure_elements = [str(e) for e in dbf.elements if str(e) not in ['VA', '/-', 'ELECTRON_GAS']]
    pure_elements.sort() # Ensures consistent ordering
    x_comp, y_comp, z_comp = pure_elements # Unpack the remaining three into your variables
    # Define the full backend components list for the solver (Must include "VA")
    solver_components = pure_elements + ['VA']
    if phases is None:
        phases = [p for p in dbf.phases.keys() if p != "NEWSIGMA"]
    # If phases are not provided, gets all possible phases programmatically with minor type error handling
    if phases is None or not isinstance(phases, list):
        phases = list(dbf.phases.keys())
        if "NEWSIGMA" in phases: phases.remove("NEWSIGMA") # Removes problematic NEWSIGMA phase from old TDB file
    # Create axes object using pycalphad built-in ternplot
    ax = ternplot(dbf, solver_components, phases, {v.T: temp, v.P: p_pa, v.X(x_comp): step_range, v.X(y_comp): step_range},
                  plot_kwargs={"tielines": True, "label_ranges": True}) # Keeps ugly "shading" in graph region
    # Add some formatting for the graph
    ax.set_title(f"Ternary Phase Diagram: {'-'.join(pure_elements)} at {temp} K", fontsize=14, pad=20)
    # Add saving capabilities
    if save:
        plt.savefig("ternary_diagram.jpg", format="jpg", dpi=300, bbox_inches="tight")
        print("Ternary phase diagram saved as ternary_diagram.jpg.")
    plt.show()

print("Functions written.")

# _________ Main Code _________
# Delete all "# RUN:" comments to show full capabilities!

### Use the retrieve_tdb function to read in the Al Li Zn database file
# RUN: db_demo = Database(retrieve_tdb("Al_Cr_Ni_Dupin_2001_TDB.TDB")) # Creates a Database object

# Use orientation function to describe the database object and data contained therein.
# RUN: orient_database(db_demo) # Calls shorter version of summary
# RUN: orient_database(db_demo, sum_length="Short") # Calls full summary

# Use binary phase diagram function; uncomment TWO lines below
# RUN: binary_components = ["AL", "NI", "VA"] # Define the components
# RUN: plot_binary_diagram(db_demo, binary_components, "AL", save=True)

# Use ternary phase diagram function; uncomment TWO lines below
# RUN: ternary_phases = ['LIQUID', 'FCC_A1', 'BCC_A2', 'SIGMA']
# RUN: plot_ternary_diagram(db_demo, temp=1500, save=True) # Current call calculates all phases at 1500

# _________ Mitch's Code: Doped Refractories for EBCs _________
# Consult: https://www.sciencedirect.com/science/article/pii/S0364591613001065?via%3Dihub
# Use this as opportunity to understand relevant .tdb file

# Retrieve the path to relevant database and then create a Database object out of it
Gd_La_Zr_db_path = retrieve_tdb("mmc1.TDB")
Gd_La_Zr_db = Database(Gd_La_Zr_db_path)

# Print Database functions and contents using orient database function
orient_database(Gd_La_Zr_db, sum_length="short")

# Create a list of all the phases in one line
Gd_La_Zr_phases = sorted(Gd_La_Zr_db.phases.keys())
# Remove metallic and intermetallic phases
Gd_La_Zr_phases = [p for p in Gd_La_Zr_phases if p not in ['BCC_A2', 'FCC_A1', 'HCP_A3', 'DHCP', 'CBCC_A12', 'CUB_A13', 'LAVES_C15', 'OMEGA', 'ORTHORHOMBIC_A20']]
# Remove gas and generic liquid (preferring IONIC_LIQ)
Gd_La_Zr_phases = [p for p in Gd_La_Zr_phases if p not in ['GAS', 'LIQUID']]
# Debug: Print the list to visualize it
# print(Gd_La_Zr_phases)

# Create a list of components to use for binary phase diagrams
Zr_La_comps = ["LA", "ZR", "O", "VA"]
Zr_Gd_comps = ["GD", "ZR", "O", "VA"]

# Debug: Stress Test
# import time
# # Test each phase individually to find the "Staller"
# for p in Gd_La_Zr_phases:
#     print(f"Testing phase: {p}...", end=" ", flush=True)
#     start = time.time()
#     try:
#         # Use calculate() to evaluate just the energy of THIS phase
#         calculate(Gd_La_Zr_db, Zr_Gd_comps, [p], T=2000, P=101325, pdens=2)
#         print(f"Passed in {time.time()-start:.2f}s")
#     except Exception as e:
#         print(f"FAILED: {e}")

# Replicate binary phase diagrams from Fig 1 and 2 from Ref. above
# Since I need to map the oxygen concentration, we cannot simply invoke binplot/function for alloys above

# x varies from 0 (ZrO2) to 0.4 (approx RE2O3 limit)
x_RE = np.linspace(0, 0.4, 100)
eq_results = []

print("Running point-by-point equilibrium with MU(O) buffer...")
for conc_RE in x_RE:
    # Use v.MU('O') instead of v.X('O') to avoid strict stoichiometry stalls
    conds_RE = {
        v.P: 101325,
        v.T: 2000,
        v.X('GD'): conc_RE,
        v.MU('O'): -200000,  # Buffer oxygen chemical potential
        v.N: 1
    }

    # Running calculation
    eq_result = equilibrium(Gd_La_Zr_db, Zr_Gd_comps, Gd_La_Zr_phases, conds_RE, calc_opts={'pdens': 20})
    eq_results.append(eq_result)

# Combine results
plot_results = xr.concat(eq_results, dim="X_GD", join="outer")
plot_results.coords['X_GD'] = x_RE
print(f"Calculated X(O) at first point: {plot_results.X.sel(component='O').values.flatten()[0]}")

# Extract and clean stable phases
stable_phases = np.unique(plot_results.Phase.values.astype(str))
print("Found phases:", stable_phases)
stable_phases = [p for p in stable_phases if p not in ['', 'nan', 'None']]

# Plotting
fig, ax = plt.subplots(figsize=(10, 6))
phase_handles, phasemap = phase_legend(stable_phases)

for phase_name in stable_phases:
    # 1. Filter and sum phase fractions (NP)
    data = plot_results.NP.where(plot_results.Phase == phase_name).sum(dim='vertex').fillna(0)

    # 2. Collapse all dimensions of size 1 (P, T, etc.)
    # This turns (1, 1, 100) into (100,)
    y_values = data.values.squeeze()

    # 3. Double check that we have a 1D array to match x_RE
    if y_values.ndim == 0:  # Handles cases where only one point exists
        y_values = [y_values]

    ax.plot(x_RE, y_values, label=phase_name, color=phasemap[phase_name], lw=2)

ax.set_title(f"Isotherm at 2000 K (MU(O) = -200 kJ/mol)")
ax.set_xlabel("Mole Fraction Gd")
ax.set_ylabel("Phase Fraction (NP)")
ax.set_ylim(0, 1.1)
ax.legend(loc='best')
plt.grid(True, alpha=0.3)
plt.show()
