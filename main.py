### CALPHAD Practice
# Mitchell Shapiro-Albert
# April 29, 2026

# _________ Imports _________
print("Beginning code...")

# Handling file imports
from pathlib import Path
# For graphing
import matplotlib.pyplot as plt
# Main package to be used
from pycalphad import Database, Workspace, binplot, ternplot, calculate
import pycalphad.variables as v # Import specific submodules as needed
from pycalphad.property_framework.metaproperties import IsolatedPhase
from pycalphad.property_framework.metaproperties import DormantPhase
from pycalphad.property_framework.tzero import T0
# To remove distracting warning messages; see note below
import warnings
# To allow for timing functions/processes
import time
import threading

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

# BETA: Write a function to plot energy curves for several phases, based on work from Anna Cassanelli
# Framework is correct, but incompatible with current TDB?
def plot_energy_curves(wrksp, x_component, save=False):
    fig = plt.figure() # Standard matplotlib
    ax = fig.add_subplot()
    x = wrksp.get(v.X(x_component)) # Retrieving the array of x_component mole fractions generated by the workspace.
    ax.set_xlabel(f"{v.X(x_component).display_name} [{v.X(x_component).display_units}]") #Labels the x‑axis
    for phase_name in wrksp.phases:  # Iterate through each phase
        # Remove copy block which caused infinite loop
        prop = IsolatedPhase(phase_name, wrksp)(f'GM({phase_name})') # Computes free energy for the phase across the composition range.
        prop.display_name = f'GM({phase_name})' # Renames the property so the legend looks clean.
        ax.plot(x, wrksp.get(prop), label=prop.display_name) # wks2.get(prop) evaluates free energy over the workspace grid
    # Add some formatting for the graph
    ax.set_title(f"Free Energy Curve for Phases: {', '.join(wrksp.phases)}", fontsize=14, pad=20)
    ax.set_ylabel("Gibbs Energy [J/mol]")
    ax.legend()
    # Add saving capabilities
    if save:
        plt.savefig("free_energy_curve.jpg", format="jpg", dpi=300, bbox_inches="tight")
        print("Free energy curve saved as free_energy_curve.jpg")
    plt.show()

# BETA: Write a function to more quickly plot energy curves for several phases
# Not debugged! Function gives LIQUID phase only with parabola.
def fast_plot_energy_curves(wrksp, x_component, save=False):
    # Extract data from the existing workspace object
    dbf = wrksp.database
    comps = wrksp.components
    phases = wrksp.phases
    # Get current temperature and pressure from the workspace conditions
    temp = wrksp.conditions.get(v.T)
    pres = wrksp.conditions.get(v.P)
    # Standard matplotlib
    fig = plt.figure()
    ax = fig.add_subplot()
    for phase_name in phases:
        # 'calculate' generates the energy surface for the phase
        # It's significantly faster than the Workspace's IsolatedPhase
        res = calculate(dbf, comps, phase_name, T=temp, P=pres, output='GM')
        # Extract the mole fraction of the specific component and the Gibbs energy
        # We use .squeeze() to remove any single-item dimensions (like T or P)
        comp_name = x_component.species.name
        x_vals = res.X.sel(component=comp_name).squeeze()
        y_vals = res.GM.squeeze()
        # Note: calculate often returns points that need sorting for a clean line plot
        sort_indices = x_vals.argsort()
        ax.plot(x_vals[sort_indices], y_vals[sort_indices], label=phase_name)
    ax.set_xlabel(f"X({x_component})")
    ax.set_ylabel("Gibbs Energy [J/mol]")
    ax.set_title(f"Fast Free Energy Curves at {temp} K", fontsize=14, pad=20)
    ax.legend()
    if save:
        plt.savefig("fast_free_energy_curve.jpg", format="jpg", dpi=300, bbox_inches="tight")
        print("Free energy curve saved as fast_free_energy_curve.jpg")
    plt.show()

print("Functions written.")

# _________ Main Code _________
# Delete all # RUN:" comments to show full capabilities!

### Use the retrieve_tdb function to read in the Al Li Zn database file
db_demo_path = retrieve_tdb("Al_Cr_Ni_Dupin_2001_TDB.TDB") # Finds the path to the database file
db_demo = Database(db_demo_path) # Creates a Database object

# Use orientation function to describe the database object and data contained therein.
orient_database(db_demo) # Calls shorter version of summary
# RUN: orient_database(db_demo, sum_length="Short") # Calls full summary

# Use binary phase diagram function; uncomment TWO lines below
# RUN: binary_components = ["AL", "NI", "VA"] # Define the components
# RUN: plot_binary_diagram(db_demo, binary_components, "AL", save=True)

# Use ternary phase diagram function; uncomment TWO lines below
# RUN: ternary_phases = ['LIQUID', 'FCC_A1', 'BCC_A2', 'SIGMA']
# RUN: plot_ternary_diagram(db_demo, temp=1500, save=True) # Current call calculates all phases at 1500

# BETA: Call energy curve functionality (based off work by Anna Cassanelli)
# Unlike above functions, this utilizes a Workspace object in pycalphad – much simpler!
# workspace_demo = Workspace(db_demo_path,
#                           ["AL", "NI"], ['LIQUID', 'FCC_A1', 'BCC_A2'],
#                           {v.X("AL"):(0,1,0.1), v.T: 1000, v.P: 101325, v.N: 1})
# plot_energy_curves(workspace_demo, v.X("AL"), save=True)
# fast_plot_energy_curves(workspace_demo, v.X("AL"), save=False) ### Not debugged

# _________ Mitch's Code: Doped Refractories for EBCs _________
# To do: Use "mmc2.TDB (Consult: https://www.sciencedirect.com/science/article/pii/S0364591613001065?via%3Dihub)
# Test orient database, binary, and ternary with this more interesting TDB
# If errors are consistent, spin ceramics_demo.py with no functionalization and add free energy capabilities
# If errors are minimal, try to resolve BETA energy curve calls