# Necessary imports
from ase.build import fcc100, add_adsorbate
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.io import write
from ase import Atoms


# Function to create a butane molecule
def create_butane():
    # Create the butane molecule
    butane = Atoms('C4H10', positions=[
        [0.0, 0.0, 0.0],  # C1
        [1.54, 0.0, 0.0],  # C2
        [2.54, 1.54, 0.0],  # C3
        [4.08, 1.54, 0.0],  # C4
        [0.0, -0.78, -0.78],  # H1
        [0.0, 0.78, -0.78],  # H2
        [-0.78, 0.0, 0.78],  # H3
        [1.54, -0.78, -0.78],  # H4
        [1.54, 0.78, -0.78],  # H5
        [2.54, 2.32, -0.78],  # H6
        [2.54, 1.54, 0.78],  # H7
        [4.08, 2.32, -0.78],  # H8
        [4.08, 1.54, 0.78],  # H9
        [4.86, 1.54, 0.0],  # H10
    ])
    butane.calc = EMT()
    return butane


# Step 1: Create the (100) silver surface (6x6x4 supercell with 10 Ångströms vacuum)
print("Erstellen der Silberoberfläche...")
surface = fcc100('Ag', size=(6, 6, 4), vacuum=10.0)
surface.calc = EMT()

# Fix the bottom two layers of the surface
mask = [atom.tag > 1 for atom in surface]
surface.set_constraint(FixAtoms(mask=mask))

# Step 2: Perform a conformer search for butane
print("Erstellen von Butan-Konformern...")
butane = create_butane()

# Optimize the butane molecule
dyn = BFGS(butane)
dyn.run(fmax=0.05)

# Generate different conformers by rotating around the C-C bonds
conformers = [butane.copy() for _ in range(3)]
angles = [0, 120, 240]  # Different angles for different conformers

for i, conf in enumerate(conformers):
    conf.rotate(angles[i], 'z')
    # Additional optimization after rotation
    conf.calc = EMT()
    dyn = BFGS(conf)
    dyn.run(fmax=0.05)

# Step 3: Place the adsorbate approximately in the middle of the supercell along the x and y axes
print("Platzierung der Konformer auf der Oberfläche...")

# Calculate the center of the surface in the x and y directions
center_x = (surface.positions[:, 0].max() + surface.positions[:, 0].min()) / 2
center_y = (surface.positions[:, 1].max() + surface.positions[:, 1].min()) / 2

# Add each conformer to the surface and save as CIF
for i, conf in enumerate(conformers):
    adsorbate = conf.copy()
    add_adsorbate(surface, adsorbate, height=2.0, position=(center_x, center_y))

    # Save the combined system as a CIF file
    filename = f'silver_butane_conformer_{i + 1}.cif'
    write(filename, surface)
    print(f"Datei '{filename}' erfolgreich erstellt.")

print("Alle CIF-Dateien wurden erfolgreich erstellt.")
