import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import pca, rms, align, diffusionmap
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.analysis.dihedrals import Dihedral
from tqdm import tqdm
from itertools import combinations
from sklearn.metrics import silhouette_score, normalized_mutual_info_score


def compute_dihedrals(universe):
    """
    计算轨迹中每个残基的 phi/psi/chi1
    不存在的 chi 用 np.nan 表示
    返回 DataFrame，列：Frame, Resid, Resname, Phi, Psi, Chi1
    """
    u = universe.copy()
    residues = list(u.residues)

    # 先收集每个残基的选择列表
    phi_list, psi_list, chi1_list = [], [], []
    resids, resnames = [], []

    for res in residues:
        phi_list.append(res.phi_selection())
        psi_list.append(res.psi_selection())
        chi1_list.append(res.chi1_selection())
        resids.append(res.resid)
        resnames.append(res.resname)

    # 只保留非 None 的选择，记录索引
    phi_indices = [i for i, sel in enumerate(phi_list) if sel is not None]
    psi_indices = [i for i, sel in enumerate(psi_list) if sel is not None]
    chi1_indices = [i for i, sel in enumerate(chi1_list) if sel is not None]

    # 创建 Dihedral 对象并计算
    phi_angles = np.full((len(u.trajectory), len(residues)), np.nan)
    psi_angles = np.full((len(u.trajectory), len(residues)), np.nan)
    chi1_angles = np.full((len(u.trajectory), len(residues)), np.nan)

    if phi_indices:
        phi_dih = Dihedral([phi_list[i] for i in phi_indices])
        phi_dih.run()
        phi_angles[:, phi_indices] = phi_dih.angles

    if psi_indices:
        psi_dih = Dihedral([psi_list[i] for i in psi_indices])
        psi_dih.run()
        psi_angles[:, psi_indices] = psi_dih.angles

    if chi1_indices:
        chi1_dih = Dihedral([chi1_list[i] for i in chi1_indices])
        chi1_dih.run()
        chi1_angles[:, chi1_indices] = chi1_dih.angles

    # 组装 DataFrame
    records = []
    for frame_idx in range(len(u.trajectory)):
        for i, res in enumerate(residues):
            records.append({
                "Frame": frame_idx,
                "Resid": resids[i],
                "Resname": resnames[i],
                "Phi": phi_angles[frame_idx, i],
                "Psi": psi_angles[frame_idx, i],
                "Chi1": chi1_angles[frame_idx, i]
            })

    df = pd.DataFrame(records)
    df["Reslabel"] = df["Resname"] + "_" + df["Resid"].astype(str)
    return df
    

def compute_RMSF(universe):
    u = universe.copy()
    align.AlignTraj(u, u, select="protein and name CA", in_memory=True).run()
    # 通常选蛋白主链或Cα原子
    ca = u.select_atoms("protein and name CA")

    # 计算 RMSF
    rmsf = rms.RMSF(ca).run()
    df_rmsf = pd.DataFrame(rmsf.rmsf)
    df_rmsf.columns = ["rmsf"]
    df_rmsf["resid"] = np.arange(1, df_rmsf.shape[0]+1, 1)
    return df_rmsf


def compute_PCA(universe):
    u = universe.copy()
    ca = u.select_atoms("name CA")
    # 2. 对齐轨迹（以去除整体运动）
    align.AlignTraj(u, u, select="name CA", in_memory=True).run()

    # 3. 运行 PCA
    pca_analysis = pca.PCA(u, select="name CA", align=False).run()

    # 4. 输出特征值（方差贡献）
    variance = pca_analysis.variance # (501, )
    # print("Explained variance (前5个主成分):")
    # print(variance[:5])

    eigvecs = pca_analysis.p_components  # shape (501, 501)

    # 5. 投影轨迹到前两个主成分
    proj = pca_analysis.transform(ca)
    df_proj = pd.DataFrame(proj)
    df_proj = df_proj.iloc[:,:3]
    df_proj.columns = ["PC1", "PC2", "PC3"]
    df_proj["frame"] = df_proj.index
    return df_proj

def compute_DCCM(universe):
    # 2. 选择主链原子（或CA原子）
    u = universe.copy()
    atom_selection = u.select_atoms("name CA")

    # 3. 对齐轨迹（去除整体平移和旋转）
    aligner = align.AlignTraj(u, u, select="name CA", in_memory=True).run()

    # 4. 计算均值位置并计算每帧的偏移 Δr
    coords = atom_selection.positions
    n_atoms = atom_selection.n_atoms
    coords_all = np.array([atom_selection.positions for ts in u.trajectory])
    mean_coords = coords_all.mean(axis=0)
    delta_r = coords_all - mean_coords

    # 5. 计算 DCCM
    n_frames = coords_all.shape[0]
    dccm = np.zeros((n_atoms, n_atoms))

    for i in range(n_atoms):
        for j in range(n_atoms):
            rij = np.sum(delta_r[:, i, :] * delta_r[:, j, :], axis=1)
            dccm[i, j] = np.mean(rij)
            
    # 归一化
    norm = np.sqrt(np.mean(np.sum(delta_r**2, axis=2), axis=0))
    dccm = dccm / np.outer(norm, norm) # shape (N, N)
    df_dccm = pd.DataFrame(dccm, index=np.arange(1, dccm.shape[0]+1, 1), columns=np.arange(1, dccm.shape[1]+1, 1))
    return df_dccm


def compute_contact_frequency(
    universe: mda.Universe,
    df_ref: pd.DataFrame,
    dist: float = 5.0,
    heavy_atoms: list = ["C","N","O","S"],
):
    """
    从 MD 轨迹构建残基接触图
    u: MDAnalysis Universe 对象
    df_ref: 参考 DataFrame，包含 Resid, Resname, Category, Reslabel 列
    dist: 接触距离阈值
    df_weight: 可选，残基间权重矩阵 DataFrame
    n_jobs: 并行计算数量（暂未实现）
    heavy_atoms: 用于计算接触的重原子名称列表
    sse_weight: SSE 权重（暂未实现）
    
    返回:
    graph_residue: NetworkX 图对象
    df_contact: 残基接触频率 DataFrame
    """
    # residues 选择部分
    u = universe.copy()
    residues = [res for res in u.residues if res.resid in df_ref['Resid'].values]
    n_res = len(residues)
    # resid_index = {res.resid: i for i, res in enumerate(residues)}

    # 初始化计数矩阵
    counts = np.zeros((n_res, n_res), dtype=int)
    # 遍历轨迹
    for ts in tqdm(u.trajectory, desc="Computing contact frequency"):
        coords_list = []
        for res in residues:
            atoms = res.atoms.select_atoms("name " + " ".join(heavy_atoms))
            coords_list.append(atoms.positions)
        
        # 计算残基对
        for i, j in combinations(range(n_res), 2):
            if len(coords_list[i]) == 0 or len(coords_list[j]) == 0:
                continue
            dists = distance_array(coords_list[i], coords_list[j], box=None, backend='serial')
            if np.any(dists <= dist):
                counts[i, j] += 1
                counts[j, i] += 1  # 对称

    n_frames = len(u.trajectory)

    # 转换为频率矩阵
    freq_matrix = counts / n_frames

    # 构建 DataFrame，行列索引为 Resid 编号或残基名
    # res_labels = [f"{res.resname}{res.resid}" for res in residues]
    df_contact = pd.DataFrame(freq_matrix,
                              index=residues, columns=residues)

    return df_contact

def compute_NMI(df, obj='Phi', n_bins=72, average_method="geometric"):
    """
    计算残基间 NMI
    df: DataFrame, 包含 Frame, Resid, angle 列
    angle: 'Phi' 或 'Psi'
    n_bins: 离散化 bin 数
    """
    residues = df['Resid'].unique()
    n_res = len(residues)
    
    def safe_digitize(data, bins):
        """
        安全版的 np.digitize:
        - 对 NaN 保持为 NaN（不分类）
        - 其他值正常分类
        """
        data = np.asarray(data)
        
        # 创建输出数组，先填 NaN
        result = np.full_like(data, fill_value=np.nan, dtype=float)
        
        # 找出非 NaN 元素
        mask = ~np.isnan(data)
        
        # 只对非 NaN 部分做 digitize
        result[mask] = np.digitize(data[mask], bins) - 1  # -1 保证从 0 开始
        
        return result
    
    # 离散化函数
    def discretize(series):
        # 将角度映射到 [0, n_bins-1]
        bins = np.linspace(-180, 180, n_bins+1)
        return safe_digitize(series, bins)
    
    
    # 构建残基 x 时间矩阵
    matrix_angle = []
    for res in residues:
        s = df[df['Resid']==res][obj].values
        matrix_angle.append(discretize(s))
    matrix_angle = np.array(matrix_angle)  # shape: (n_res, n_frames)
    
    # 计算 NMI 矩阵
    nmi_matrix = np.zeros((n_res, n_res))
    for i in range(n_res):
        x = matrix_angle[i,:]
        mask_x = ~(np.isnan(x))
        x_valid = x[mask_x]
        if len(x_valid) == 0:
            nmi_matrix[i, :] = np.nan
            nmi_matrix[:, i] = np.nan
            continue
        
        for j in range(i, n_res):
            y = matrix_angle[j,:]
            # print("x: ", x)
            # print("y: ", y)
            mask = ~(np.isnan(x) | np.isnan(y))
            x_valid, y_valid = x[mask], y[mask]
            if len(x_valid) == 0 or len(y_valid) == 0:
                nmi_matrix[i,j] = np.nan
                nmi_matrix[j,i] = np.nan
                continue
            
            nmi = normalized_mutual_info_score(x_valid, y_valid, average_method=average_method)
            nmi_matrix[i,j] = nmi
            nmi_matrix[j,i] = nmi  # 对称
    
    df_nmi = pd.DataFrame(nmi_matrix, index=residues, columns=residues)
    return df_nmi

def compute_weighted_average(dfs, weights=None, how='inner'):
    """
    对多个 DataFrame 进行加权平均，自动忽略 NaN。
    
    参数:
      dfs: list[pd.DataFrame]
          要合并的 DataFrame 列表
      weights: list[float] 或 None
          每个 DataFrame 的权重；若为 None，则均分
      how: str
          对齐方式，'inner'（交集）或 'outer'（并集）
    
    返回:
      pandas.DataFrame
    """
    if not dfs:
        raise ValueError("dfs 不能为空")
    
    n = len(dfs)
    if weights is None:
        weights = [1.0] * n
    if len(weights) != n:
        raise ValueError("weights 长度必须与 dfs 数量一致")
    
    # 对齐所有 DataFrame
    df_aligned = dfs[0].copy()
    for df in dfs[1:]:
        df_aligned, df = df_aligned.align(df, join=how)
    
    # 将所有 DataFrame 对齐到相同的索引/列
    aligned_dfs = [df.align(df_aligned, join='outer')[0] for df in dfs]

    # 转换为三维数组，方便向量化计算
    stacked = np.stack([df.to_numpy() for df in aligned_dfs], axis=-1)
    weights = np.array(weights).reshape(1, 1, -1)

    # 计算加权平均（忽略 NaN）
    valid_mask = ~np.isnan(stacked)
    weighted_sum = np.nansum(stacked * weights, axis=-1)
    weight_sum = np.nansum(weights * valid_mask, axis=-1)
    result = weighted_sum / np.where(weight_sum == 0, np.nan, weight_sum)

    # 转回 DataFrame
    return pd.DataFrame(result, index=df_aligned.index, columns=df_aligned.columns)


def compute_edge_weight(df_NMI, df_contact):
    """
    整合 DCCM 和 NMI 矩阵，生成新的相关性矩阵
    dccm: np.ndarray, DCCM 矩阵
    df_NMI: pd.DataFrame, NMI 矩阵
    返回: pd.DataFrame, 整合后的相关性矩阵
    """
    df_weight = df_NMI.values * df_contact.values
    df_weight = pd.DataFrame(df_weight, index=df_NMI.index, columns=df_NMI.columns)
    return df_weight
