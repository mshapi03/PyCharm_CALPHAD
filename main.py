### CALPHAD Practice
# Mitchell Shapiro-Albert
# April 23, 2026

# _________ Imports _________
print("Beginning code...")

from pathlib import Path
import matplotlib.pyplot as plt
from pycalphad import Database, binplot
import pycalphad.variables as v # Import specific submodules as needed
import warnings

print("Imports successful.")

# _________ Warning Handling _________
# NOTE: you should remove the warning handling if you are using this for research purposes! These suppress warning
# messages that you should be aware of as a user. E.g., if you are calculating molar volume or magnetic heat capacity,
# these warnings will tell you that the .tdb file you are using is not a good candidate.
warnings.filterwarnings("ignore", message=".*Type definitions using IF/THEN logic.*")
warnings.filterwarnings("ignore", message=".*no corresponding TYPE_DEFINITION line was found.*")

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

# Write a function to plot a binary phase diagram given a retrieved_tdb object, components, and a temperature range;
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

print("Functions written.")

# _________ Main Code _________

# Use the retrieve_tdb function to read in the Al Li Zn database file
db_current = Database(retrieve_tdb("Al_Cr_Ni_Dupin_2001_TDB.TDB")) # Creates a Database object

# Use orientation function to describe the database object and data contained therein.
# orient_database(db_current, sum_length != "Full") # Calls shorter version of summary
orient_database(db_current) # Calls full summary

# Use binary phase diagram function; uncomment TWO lines below
# binary_components = ["AL", "NI", "VA"] # Define the components
# plot_binary_diagram(db_current, binary_components, "AL", save=True )

# Use ternary phase diagram function; uncomment _ lines below

