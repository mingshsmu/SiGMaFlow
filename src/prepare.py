from Bio.PDB import PDBParser
import pandas as pd
import numpy as np
import mdtraj as md
import MDAnalysis as mda
from itertools import combinations





aa_classification = {
    # ===== Acidic =====
    "ASP": "Acidic", "GLU": "Acidic",
    # 质子化形式
    "ASH": "Acidic", "GLH": "Acidic",

    # ===== Basic =====
    "LYS": "Basic", "ARG": "Basic", "HIS": "Basic",
    # Histidine tautomers
    "HID": "Basic", "HIE": "Basic", "HIP": "Basic",

    # ===== Aromatic =====
    "PHE": "Aromatic", "TYR": "Aromatic", "TRP": "Aromatic",

    # ===== Polar (uncharged) =====
    "SER": "Polar", "THR": "Polar", "ASN": "Polar", "GLN": "Polar",

    # ===== Hydrophobic =====
    "ALA": "Hydrophobic", "VAL": "Hydrophobic", "ILE": "Hydrophobic",
    "LEU": "Hydrophobic", "MET": "Hydrophobic", "PRO": "Hydrophobic",

    # ===== Special =====
    "GLY": "Special", "CYS": "Special",

    # ===== Seleno and Modified =====
    "MSE": "Hydrophobic",   # Selenomethionine
    "CSE": "Special",       # Selenocysteine

    # ===== Phosphorylated variants =====
    "SEP": "Polar",  # Phosphoserine
    "TPO": "Polar",  # Phosphothreonine
    "PTR": "Aromatic"  # Phosphotyrosine
}


def classify_residue(resname: str, aa_dict: dict) -> str:
    """
    根据残基名分类:
    1. 直接查字典
    2. 如果不存在，尝试前三个字母
    3. 否则返回 Unknown
    """
    resname = resname.upper().strip()

    # 直接查
    if resname in aa_dict:
        return aa_dict[resname]

    # 尝试前三个字母 (避免过短字符串报错)
    if len(resname) >= 3 and resname[:3] in aa_dict:
        return aa_dict[resname[:3]]

    # fallback
    return "Unknown"


def parse_pdb(input_path):
    # 解析 PDB
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_path)
    records = []
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                resid = residue.id[1]   # 残基序号
                if "CA" in residue:     # 只取Cα
                    ca = residue["CA"]
                    x, y, z = ca.coord
                    aa_class = classify_residue(resname, aa_classification)

                    records.append({
                        "Chain": chain.id,
                        "Resid": resid,
                        "Resname": resname,
                        "Category": aa_class,
                        "CA_x": x,
                        "CA_y": y,
                        "CA_z": z
                    })

    df = pd.DataFrame(records)
    df["Reslabel"] = df["Resname"] + "_" + df["Resid"].astype(str)
    
    traj = md.load(input_path) 
    dssp = md.compute_dssp(traj)[0]
    pre_sse = None
    segments = []
    segment_id = 0
    for sse in dssp:
        if sse != pre_sse:
            segment_id += 1
            pre_sse = sse
        
        segments.append(segment_id)

    df['SSE'] = dssp
    df['Segment'] = segments
    
    return df


def calculate_AA_dist(df_ref_pdb):
    coord_cols=['CA_x', 'CA_y', 'CA_z']
    resid = df_ref_pdb['Resid'].values
    coords = df_ref_pdb[coord_cols].values  # shape = (N,3)

    # 生成所有残基对
    pairs = list(combinations(range(len(resid)), 2))

    # 计算距离
    distances = []
    for i,j in pairs:
        dist = np.linalg.norm(coords[i] - coords[j])
        distances.append(dist)

    # 构建距离表
    df_AApair_dist = pd.DataFrame({
        'Resid1': [resid[i] for i,j in pairs],
        'Resid2': [resid[j] for i,j in pairs],
        'distance': distances
    })
    return df_AApair_dist

