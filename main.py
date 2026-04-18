### CALPHAD Practice
# Mitchell Shapiro-Albert
# April 16, 2026

# _________ Imports _________
print("Beginning code...")

from pathlib import Path
import matplotlib.pyplot as plt
from pycalphad import Database, binplot
import pycalphad.variables as v # Import specific submodules as needed

print("Imports successful.")

# _________ Functions _________
print("Writing functions...")

# Write a function to retrieve the .tbd file from the subdirectory
def retrieve_tdb(file_name):
    script_dir = Path(__file__).resolve()
    target_path = script_dir.parent / "TDB_Files" / f"{file_name}"
    if target_path.exists():
        return target_path
    else:
        raise FileNotFoundError(f"{target_path} does not exist")

# Write a function to fully summarize a Database object
def orient_database(dbf, type="Full"):
    print(f"{'=' * 20} IMPORTED DATABASE SUMMARY {'=' * 20}")

    # 1. Basic Components and Species
    print(f"\n[COMPONENTS & SPECIES]")
    print(f"Elements: {sorted(dbf.elements)}")
    print(f"Species:  {[s.name for s in dbf.species]}")

    if type == "Full":
        # 2. Parameters and Functions (Symbols)
        print(f"\n[SYMBOLS/FUNCTIONS]")
        functions = [sym for sym in dbf.symbols.keys()]
        print(f"Total Symbols: {len(functions)}")
        if functions:
            print(f"Sample Symbols: {functions[:5]}...")

    # 3. Detailed Phase Information
    print(f"\n[PHASE DESCRIPTIONS]")
    phase_names = sorted(dbf.phases.keys())
    print(f"Total Phases: {len(phase_names)}")

    for name in phase_names:
        phase = dbf.phases[name]
        # Extract sublattice and constituent information
        sublattices = phase.sublattices
        # Constituents are stored as a tuple of species for each sublattice
        constituents = [[s.name for s in c] for c in phase.constituents]

        print(f"\nPhase: {name}")
        if type == "Full":
            print(f"  - Sublattices: {sublattices}")
            for i, (sites, species) in enumerate(zip(sublattices, constituents)):
                print(f"    Sublattice {i + 1} ({sites} sites): {', '.join(species)}")

            # Check for model hints (e.g., if it's a magnetic or ordered phase)
            if hasattr(phase, 'model_hints') and phase.model_hints:
                print(f"    Hints: {list(phase.model_hints.keys())}")

    print(f"\n{'=' * 60}")

# Write a function to plot a binary phase diagram given a file, components, and a temperature range;
# Formerly broken, currently non-existent
def plot_binary_diagram(dbf, components, temp_range=(300, 1000, 10)):
    pass

print("Functions written.")

# _________ Main Codebase _________

# Use the retrieve_tdb function to read in the Al Li Zn database file
db_current = Database(retrieve_tdb("Al_Cr_Ni_Dupin_2001_TDB.TDB")) # Creates a Database object

# Use orientation function to describe the database object and data contained therein.
# orient_database(db_current, type="Short") # Calls shorter version of summary
orient_database(db_current) # Calls full summary

# Manual Plot generation
print("Graphing...")
components = ["AL", "NI", "VA"] # Define the components
phases = db_current.phases.keys() # Get all possible phases programmatically
binplot(db_current, components, phases, {v.N: 1, v.P: 101325, v.T: (300, 2000, 10), v.X("AL"):(0, 1, 0.05)}, plot_kwargs={"tielines":False})
plt.show()

# Note: works! Functionalize this process and add a warning handling for the UserWarning that pops up.
# Then do the same with ternary diagrams :)