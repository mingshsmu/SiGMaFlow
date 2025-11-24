from pymol import cmd, cgo
cmd.show("cartoon", "polymer.protein")
cmd.show("sticks", "organic")
cmd.util.cnc("all")
cmd.select('stable_AA', 'resi 3+4+5+6+7+51+52+53+54+55+74+75+76+77+78+79+80+81+93+94+95+96+107+108+109+110+111+112+113+114+115+137+138+139+140+141+157+158+161')
cmd.select('fluc_AA', 'resi 25+26+27+28+29+30+33+34+35+59+60+61+62+63+64+65+68+147+148+167')
cmd.super('PC1_0 and stable_AA and name CA', 'PC1_2 and stable_AA and name CA')
motion_cgo = []

coord1 = cmd.get_atom_coords('PC1_0 and resi 25 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 25 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 26 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 26 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 27 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 27 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 28 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 28 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 29 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 29 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 30 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 30 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 33 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 33 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 34 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 34 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 35 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 35 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 59 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 59 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 60 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 60 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 61 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 61 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 62 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 62 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 63 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 63 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 64 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 64 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 65 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 65 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 68 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 68 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 147 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 147 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 148 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 148 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])

coord1 = cmd.get_atom_coords('PC1_0 and resi 167 and name CA')
coord2 = cmd.get_atom_coords('PC1_2 and resi 167 and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, 0.5, 0.05, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 0.9921568627450981, 0.9058823529411765, 0.1450980392156863, 1, 1])
cmd.load_cgo(motion_cgo, 'motion')
cmd.zoom('polymer.protein')
