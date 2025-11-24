from chimerax.core.commands.atomspec import AtomSpecArg
from chimerax.core.commands import run

topRadius, radius = (0.05, 0.5)
color = '#fde725'

run(session, 'align #2/X:3,4,5,6,7,51,52,53,54,55,74,75,76,77,78,79,80,81,93,94,95,96,107,108,109,110,111,112,113,114,115,137,138,139,140,141,157,158,161@CA toAtoms #1/X:3,4,5,6,7,51,52,53,54,55,74,75,76,77,78,79,80,81,93,94,95,96,107,108,109,110,111,112,113,114,115,137,138,139,140,141,157,158,161@CA')

atom_spec = '#1/X:25,26,27,28,29,30,33,34,35,59,60,61,62,63,64,65,68,147,148,167@CA'
parsed = AtomSpecArg.parse(atom_spec, session)
atom_spec_obj = parsed[0].evaluate(session)
atoms = getattr(atom_spec_obj, "atoms", None)
coords = [tuple(atom.scene_coord) for atom in atoms]
coords_1 = coords

atom_spec = '#2/X:25,26,27,28,29,30,33,34,35,59,60,61,62,63,64,65,68,147,148,167@CA'
parsed = AtomSpecArg.parse(atom_spec, session)
atom_spec_obj = parsed[0].evaluate(session)
atoms = getattr(atom_spec_obj, "atoms", None)
coords = [tuple(atom.scene_coord) for atom in atoms]
coords_2 = coords

for coord1, coord2 in zip(coords_1, coords_2):
    cmd = (f'shape cone fromPoint {coord1[0]}, {coord1[1]}, {coord1[2]} '
           f'toPoint {coord2[0]}, {coord2[1]}, {coord2[2]} '
           f'radius {radius} topRadius {topRadius} color {color} name motion')
    run(session, cmd)
