# Ligand-centred pocket visualisation in PyMOL
# Author: Dr. Jineetkumar Gawad
#
# Run in a fresh PyMOL session.
# Set PyMOL's working directory to the folder containing your inputs.
# Edit the two filenames below to match your files.
# Use a single receptor structure and a single ligand pose.
#
# The surface represents a local residue environment.
# Distance-based neighbours are not necessarily direct interactions.

# 1. Load inputs
load protein.pdb, receptor
load ligand.sdf, ligand

# 2. Select complete residues within 5 Angstrom of ligand heavy atoms
# Change 5.0 below if a different distance cutoff is required.
select pocket_residues, byres ((receptor and polymer.protein and not hydro) within 5.0 of (ligand and not hydro))

# 3. Create a separate object for the local surface
create pocket_surface, pocket_residues
hide everything, receptor
hide everything, ligand
hide everything, pocket_surface

# 4. Display the grey local surface
show surface, pocket_surface
set surface_color, grey70, pocket_surface
set transparency, 0.50, pocket_surface

# 5. Display ligand and neighbouring residues
show sticks, ligand and not hydro
color yellow, ligand
show sticks, pocket_residues and not hydro
color green, pocket_residues
set stick_radius, 0.15, receptor
set stick_radius, 0.22, ligand

# 6. Label residues at their alpha-carbon atoms
label pocket_residues and name CA, "%s:%s%s" % (chain, resn, resi)
set label_color, black
set label_size, 14
set label_font_id, 7

# 7. Set the view
bg_color white
orient ligand
zoom ligand, 7
deselect

# 8. Print and export unique residue identifiers
# Uses all selected atoms, so residues missing a CA are still included.
python
from pymol import cmd
import csv

residues = set()
cmd.iterate(
    "pocket_residues",
    "residues.add((model, segi, chain, resn, resi))",
    space={"residues": residues}
)

with open("pocket_residues.csv", "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["object", "segment", "chain", "residue_name", "residue_id"])
    for residue in sorted(residues):
        writer.writerow(residue)
        print(residue)

print("Selected residues:", len(residues))
print("Saved residue list: pocket_residues.csv")
python end

# 9. Optional exports
# After adjusting the view, remove the leading # to enable these commands.
# Existing files with these names will be overwritten.
# png pocket_view.png, width=2400, height=1800, dpi=300, ray=1
# save pocket_session.pse
