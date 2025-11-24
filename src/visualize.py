import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import json

# plt.rcParams.update({
#     "font.family": "Arial",
#     "text.usetex": False,
#     "pdf.fonttype": 42,
#     "ps.fonttype": 42,
#     "font.size": 20,
#     "axes.titlesize": 20,
#     "axes.labelsize": 20,
#     "xtick.labelsize": 20,
#     "ytick.labelsize": 20,
#     "legend.fontsize": 20,
#     "legend.title_fontsize": 20
# })


def visualize_RMSF_pymol(df_rmsf, outfile, cmap='coolwarm',midpoint=None, half_range=None):
    df = df_rmsf.copy()
    if midpoint is None:
        midpoint = df["rmsf"].median()
    if half_range is None:
        half_range = df["rmsf"].std() / 2
        
    fold = float(50 / half_range)
    limit_lower = -50
    limit_upper = 50
    cmap_colors = plt.get_cmap(cmap)
    df["color_index"] = np.round((df["rmsf"] - midpoint) * fold)
    df["color_index"] = df["color_index"].clip(limit_lower, limit_upper)
    df["color_index"] = df["color_index"] + 50
    df["color_index"] = df["color_index"].astype(int)
    
    def get_rgb(idx):
        rgb = cmap_colors(idx/100)[:3]
        rgb = f"[{rgb[0]:.3f}, {rgb[1]:.3f}, {rgb[2]:.3f}]"
        return rgb

    df["color_rgb"] = df["color_index"].apply(lambda x: pd.Series(get_rgb(x)))
    with open(outfile, "w") as f:
        for index, row in df.iterrows():
            f.write(f"set_color color_{row['color_index']}, {row['color_rgb']}\n")
            f.write(f"color color_{row['color_index']}, resi {row['resid']}\n")

def visualize_RMSF_bar(outfile=None, cmap='coolwarm', midpoint=None, half_range=None, fig_size=(0.5, 8), n_break=7):
    cmap = plt.get_cmap(cmap)
    vmin = midpoint - half_range
    vmax = midpoint + half_range
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    fig, ax = plt.subplots(figsize=fig_size)
    fig.subplots_adjust(bottom=0.5)
    cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation="vertical", ticks=np.linspace(vmin, vmax, n_break))
    cb.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
    if outfile is not None:
        plt.savefig(f"{outfile}.png", dpi=600, bbox_inches="tight")
        plt.savefig(f"{outfile}.pdf", bbox_inches="tight")
    plt.show()

def visualize_RMSF_chimeraX(df_rmsf, outfile, cmap='coolwarm',midpoint=None, half_range=None):
    df = df_rmsf.copy()
    if midpoint is None:
        midpoint = df["rmsf"].median()
    if half_range is None:
        half_range = df["rmsf"].std() / 2
        
    fold = float(50 / half_range)
    limit_lower = -50
    limit_upper = 50
    cmap_colors = plt.get_cmap(cmap)
    df["color_index"] = np.round((df["rmsf"] - midpoint) * fold)
    df["color_index"] = df["color_index"].clip(limit_lower, limit_upper)
    df["color_index"] = df["color_index"] + 50
    df["color_index"] = df["color_index"].astype(int)
    
    def get_hex(idx):
        hex_color = cmap_colors(idx/100)
        hex_color = mcolors.to_hex(hex_color)
        return hex_color

    df["color_hex"] = df["color_index"].apply(lambda x: pd.Series(get_hex(x)))
    with open(outfile, "w") as f:
        for index, row in df.iterrows():
            f.write(f"color :{row['resid']} {row['color_hex']}\n")


def visualize_DccmCluster_pymol(dict_dccm_clusters, outfile, colors=None):
    if colors is None:
        colors = [
        [0.894, 0.102, 0.110],  # red
        [0.216, 0.494, 0.722],  # blue
        [0.302, 0.686, 0.290],  # green
        [0.596, 0.306, 0.639],  # purple
        [1.000, 0.498, 0.000],  # orange
        [1.000, 1.000, 0.200],  # yellow
        [0.651, 0.337, 0.157],  # brown
        [0.969, 0.506, 0.749],  # pink
        [0.600, 0.600, 0.600],  # gray
        [0.121, 0.470, 0.705],  # light blue
        [0.682, 0.780, 0.909],  # sky blue
        [0.200, 0.627, 0.173],  # emerald
        [0.984, 0.604, 0.600],  # salmon
        [0.890, 0.466, 0.760],  # violet
        [0.498, 0.788, 0.498],  # mint
        [0.737, 0.741, 0.133],  # olive
        [0.090, 0.745, 0.811],  # turquoise
        [0.992, 0.682, 0.380],  # peach
        [0.749, 0.357, 0.090],  # chestnut
        [0.416, 0.239, 0.604]   # indigo
        ]
    
    n_clusters = len(dict_dccm_clusters.keys())
    with open(outfile, "w") as f:
        color_names = []
        for i, rgb in enumerate(colors):
            color_name = f"cluster_col_{i+1}"
            f.write(f"set_color {color_name}, [{rgb[0]:.3f}, {rgb[1]:.3f}, {rgb[2]:.3f}]\n")
            color_names.append(color_name)
        for cluster_id, residues in dict_dccm_clusters.items():
            res_selection = "resi " + "+".join([str(r) for r in residues])
            color = color_names[cluster_id % len(color_names) - 1]
            
            f.write(f"select cluster_{cluster_id}, ({res_selection})\n")
            f.write(f"color {color}, cluster_{cluster_id}\n")
        f.write("show cartoon, all\n")
        f.write("util.cnc\n")

def visualize_DccmCluster_chimeraX(dict_dccm_clusters, outfile, colors=None):
    if colors is None:
        colors = [
        [0.894, 0.102, 0.110],  # red
        [0.216, 0.494, 0.722],  # blue
        [0.302, 0.686, 0.290],  # green
        [0.596, 0.306, 0.639],  # purple
        [1.000, 0.498, 0.000],  # orange
        [1.000, 1.000, 0.200],  # yellow
        [0.651, 0.337, 0.157],  # brown
        [0.969, 0.506, 0.749],  # pink
        [0.600, 0.600, 0.600],  # gray
        [0.121, 0.470, 0.705],  # light blue
        [0.682, 0.780, 0.909],  # sky blue
        [0.200, 0.627, 0.173],  # emerald
        [0.984, 0.604, 0.600],  # salmon
        [0.890, 0.466, 0.760],  # violet
        [0.498, 0.788, 0.498],  # mint
        [0.737, 0.741, 0.133],  # olive
        [0.090, 0.745, 0.811],  # turquoise
        [0.992, 0.682, 0.380],  # peach
        [0.749, 0.357, 0.090],  # chestnut
        [0.416, 0.239, 0.604]   # indigo
        ]
    
    n_clusters = len(dict_dccm_clusters.keys())
    color_hex = [mcolors.to_hex(color) for color in colors]
    with open(outfile, "w") as f:
        for cluster_id, residues in dict_dccm_clusters.items():
            res_selection = ":" + ",".join([str(r) for r in residues])
            color = color_hex[cluster_id % len(color_hex) - 1]
            f.write(f"color {res_selection} {color}\n")

def visualize_path_preprocess(df_paths, df_edges, df_ref_pdb, outfile=None, cluster_by="path_id", colors=None, radius_range=(0.1, 0.5)):
    df_paths = df_paths.copy()
    df_ref_pdb = df_ref_pdb.copy()
    df_edges = df_edges.copy()
    df_paths["path_id"] = "path_" + df_paths["Resid1"].astype(str) + "_" + df_paths["Resid2"].astype(str)
    if cluster_by == "cluster":
        df_paths["cluster"] = "cluster_" + df_paths["cluster"].astype(str) + "_" + df_paths["path_id"]
    
    df = df_ref_pdb[["Resid", "Resname", "CA_x", "CA_y", "CA_z"]]
    df_edges = df_edges.merge(
        df.rename(columns={
            "Resid": "Resid1",
            "Resname": "Resname1",
            "CA_x": "CA_x1",
            "CA_y": "CA_y1",
            "CA_z": "CA_z1",
        }),
        on="Resid1", how="left"
    )
    df_edges = df_edges.merge(
        df.rename(columns={
            "Resid": "Resid2",
            "Resname": "Resname2",
            "CA_x": "CA_x2",
            "CA_y": "CA_y2",
            "CA_z": "CA_z2",
        }),
        on="Resid2", how="left"
    )
    
    path_edges = []
    for _, row in df_paths.iterrows():
        path = row["path"]
        cluster = row[cluster_by]
        for i in range(len(path)-1):
            path_edges.append((path[i], path[i+1], cluster))
            path_edges.append((path[i+1], path[i], cluster))
    
    df_paths_edges = pd.DataFrame(path_edges, columns=["Resid1", "Resid2", "cluster"]).drop_duplicates()
    df_edges = df_edges.merge(df_paths_edges, on=["Resid1", "Resid2"], how="left")
    df_edges["coord1"] = list(zip(df_edges["CA_x1"], df_edges["CA_y1"], df_edges["CA_z1"]))
    df_edges["coord2"] = list(zip(df_edges["CA_x2"], df_edges["CA_y2"], df_edges["CA_z2"]))
    df_edges = df_edges.dropna(subset=["cluster"])
    df_edges = df_edges.drop(columns=["CA_x1", "CA_y1", "CA_z1", "CA_x2", "CA_y2", "CA_z2"])
    
    weights = df_edges["weight"].values
    weights_norm = (weights - np.min(weights)) / (np.max(weights) - np.min(weights))
    radius_min, radius_max = radius_range
    radius = radius_min + weights_norm * (radius_max - radius_min)
    df_edges["radius"] = radius
    
    
    dict_coords = {}
    list_clusters = list(set(df_edges["cluster"].tolist()))
    for i, row in df_edges.iterrows():
        cluster = row["cluster"]
        resid1 = row["Resid1"]
        resid2 = row["Resid2"]
        coord1 = row["coord1"]
        coord2 = row["coord2"]
        weight = row.get("weight", 0.0)
        radius = row.get("radius", 0.0)
        
        if cluster not in dict_coords:
            dict_coords[cluster] = []

        dict_coords[cluster].append({
            "path_idx": list_clusters.index(cluster),
            "cluster": cluster,
            "resid1": resid1,
            "resid2": resid2,
            "coord1": coord1,
            "coord2": coord2,
            "weight": weight,
            "radius": radius
        })
    
    if outfile is not None:
        with open(outfile, "w") as f:
            json.dump(dict_coords, f, indent=2)
    
    return dict_coords
    
def visualize_path_pymol(dict_coords, outfile=None, colors=None):
    if outfile is None:
        outfile = "vis_path.py"
    
    if colors is None:
        colors = [
            [0.902, 0.294, 0.208],
            [0.302, 0.733, 0.835],
            [0.000, 0.627, 0.529],
            [0.235, 0.329, 0.533],
            [0.953, 0.608, 0.498],
            [0.518, 0.569, 0.706],
            [0.569, 0.820, 0.760],
            [0.863, 0.000, 0.000],
            [0.494, 0.380, 0.282],
            [0.690, 0.612, 0.522],
            [1.000, 0.498, 0.055],
            [0.173, 0.627, 0.173],
            [0.122, 0.467, 0.706],
            [0.580, 0.404, 0.741],
            [0.549, 0.337, 0.294],
            [0.090, 0.745, 0.812],
            [0.839, 0.153, 0.157],
            [0.737, 0.741, 0.133],
            [0.168, 0.482, 0.729],
            [0.259, 0.710, 0.251],
            ]
    
    cmd_pymol = """from pymol import cmd, cgo
cmd.show("cartoon", "polymer.protein")
cmd.show("sticks", "organic")
cmd.util.cnc("all")
"""

    for cluster_id, rows in dict_coords.items():
        cluster_cgo = []
        for row in rows:
            pos1 = row["coord1"]
            pos2 = row["coord2"]
            radius = row["radius"]
            cluster_color = colors[int(row["path_idx"]) % len(colors)]
            # 注意：这里直接构造 list，不是字符串！
            cluster_cgo.append(
                [ "cgo.CYLINDER", *pos1, *pos2, radius, *cluster_color, *cluster_color ]
            )

        # 写成 Python 有效语句字符串
        cmd_pymol += f"{cluster_id} = [\n"
        for c in cluster_cgo:
            cmd_pymol += "    " + ", ".join(str(v) for v in c) + ",\n"
        cmd_pymol += "]\n"
        cmd_pymol += f"cmd.load_cgo({cluster_id}, '{cluster_id}')\n"

    cmd_pymol += "cmd.zoom('polymer.protein')\n"

    with open(outfile, "w") as f:
        f.write(cmd_pymol)
        
def visualize_path_chimeraX(dict_coords, outfile=None, colors=None):
    if outfile is None:
        outfile = "vis_path.cxc"
    
    if colors is None:
        colors = [
            [0.902, 0.294, 0.208],
            [0.302, 0.733, 0.835],
            [0.000, 0.627, 0.529],
            [0.235, 0.329, 0.533],
            [0.953, 0.608, 0.498],
            [0.518, 0.569, 0.706],
            [0.569, 0.820, 0.760],
            [0.863, 0.000, 0.000],
            [0.494, 0.380, 0.282],
            [0.690, 0.612, 0.522],
            [1.000, 0.498, 0.055],
            [0.173, 0.627, 0.173],
            [0.122, 0.467, 0.706],
            [0.580, 0.404, 0.741],
            [0.549, 0.337, 0.294],
            [0.090, 0.745, 0.812],
            [0.839, 0.153, 0.157],
            [0.737, 0.741, 0.133],
            [0.168, 0.482, 0.729],
            [0.259, 0.710, 0.251],
            ]
    
    colors = [mcolors.to_hex(c) for c in colors]
    
    cmd_chimeraX = ""

    for cluster_id, rows in dict_coords.items():
        cluster_cgo = []
        for row in rows:
            pos1 = row["coord1"]
            pos1 = ",".join(str(x) for x in pos1)
            pos2 = row["coord2"]
            pos2 = ",".join(str(x) for x in pos2)
            radius = row["radius"]
            cluster_color = colors[int(row["path_idx"]) % len(colors)]
            cmd_chimeraX += f"shape cylinder fromPoint {pos1} toPoint {pos2} radius {radius} color {cluster_color} name {cluster_id}\n"


    with open(outfile, "w") as f:
        f.write(cmd_chimeraX)        

def visualize_motion_pymol(obj1, obj2, df_rmsf, outfile, low_threshold=1, high_threshold=3,
                           cone_color="orange", cone_sep=2, cone_radius_range=(0.05, 0.5)):
    df = df_rmsf.copy()
    list_stable_AA = df.loc[df["rmsf"] < float(low_threshold)]["resid"].tolist()
    list_fluc_AA = df.loc[df["rmsf"] > float(high_threshold)]["resid"].tolist()
    list_str_stable_AA = [str(r) for r in list_stable_AA]
    list_str_fluc_AA = [str(r) for r in list_fluc_AA]
    select_stable_AA = "+".join(list_str_stable_AA)
    select_fluc_AA = "+".join(list_str_fluc_AA)
    list_fluc_AA = list_fluc_AA[::int(cone_sep)]
    
    cone_color = mcolors.to_rgb(cone_color)
    cone_color = ", ".join(str(r) for r in cone_color)
    
    cmd_pymol = """from pymol import cmd, cgo
cmd.show("cartoon", "polymer.protein")
cmd.show("sticks", "organic")
cmd.util.cnc("all")
"""
    cmd_pymol += f"cmd.select('stable_AA', 'resi {select_stable_AA}')\n"
    cmd_pymol += f"cmd.select('fluc_AA', 'resi {select_fluc_AA}')\n"
    cmd_pymol += f"cmd.super('{obj1} and stable_AA and name CA', '{obj2} and stable_AA and name CA')\n"
    cmd_pymol += "motion_cgo = []\n"
    
    for res_id in list_fluc_AA:
        cmd_pymol += f"""
coord1 = cmd.get_atom_coords('{obj1} and resi {res_id} and name CA')
coord2 = cmd.get_atom_coords('{obj2} and resi {res_id} and name CA')
motion_cgo.extend([cgo.CONE, *coord1, *coord2, {cone_radius_range[1]}, {cone_radius_range[0]}, {cone_color}, {cone_color}, 1, 1])
"""
    cmd_pymol += "cmd.load_cgo(motion_cgo, 'motion')\n"
    cmd_pymol += "cmd.zoom('polymer.protein')\n"
    
    with open(outfile, "w") as f:
        f.write(cmd_pymol)
        

def visualize_motion_chimeraX(obj1, obj2, df_rmsf, outfile, low_threshold=1, high_threshold=3,
                           cone_color="orange", cone_sep=2, cone_radius_range=(0.05, 0.5)):
    df = df_rmsf.copy()
    list_stable_AA = df.loc[df["rmsf"] < float(low_threshold)]["resid"].tolist()
    list_str_stable_AA = [str(r) for r in list_stable_AA]
    select_stable_AA = ",".join(list_str_stable_AA)
    
    list_fluc_AA = df.loc[df["rmsf"] > float(high_threshold)]["resid"].tolist()
    list_fluc_AA = list_fluc_AA[::int(cone_sep)]
    list_str_fluc_AA = [str(r) for r in list_fluc_AA]
    select_fluc_AA = ",".join(list_str_fluc_AA)
    
    cone_color = mcolors.to_hex(cone_color)
    
    cmd_chimeraX = f"""from chimerax.core.commands.atomspec import AtomSpecArg
from chimerax.core.commands import run

topRadius, radius = {cone_radius_range}
color = '{cone_color}'

run(session, 'align #2/X:{select_stable_AA}@CA toAtoms #1/X:{select_stable_AA}@CA')

atom_spec = '{obj1}/X:{select_fluc_AA}@CA'
parsed = AtomSpecArg.parse(atom_spec, session)
atom_spec_obj = parsed[0].evaluate(session)
atoms = getattr(atom_spec_obj, "atoms", None)
coords = [tuple(atom.scene_coord) for atom in atoms]
coords_1 = coords

atom_spec = '{obj2}/X:{select_fluc_AA}@CA'
parsed = AtomSpecArg.parse(atom_spec, session)
atom_spec_obj = parsed[0].evaluate(session)
atoms = getattr(atom_spec_obj, "atoms", None)
coords = [tuple(atom.scene_coord) for atom in atoms]
coords_2 = coords
"""
    cmd_chimeraX += """
for coord1, coord2 in zip(coords_1, coords_2):
    cmd = (f'shape cone fromPoint {coord1[0]}, {coord1[1]}, {coord1[2]} '
           f'toPoint {coord2[0]}, {coord2[1]}, {coord2[2]} '
           f'radius {radius} topRadius {topRadius} color {color} name motion')
    run(session, cmd)
"""

    with open(outfile, "w") as f:
        f.write(cmd_chimeraX)