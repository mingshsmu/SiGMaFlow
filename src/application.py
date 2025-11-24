from collections.abc import Iterable
import itertools
import numpy as np
import pandas as pd
import networkx as nx
from sklearn.cluster import SpectralClustering
import MDAnalysis as mda
from MDAnalysis.analysis import pca, rms, align, diffusionmap
import os

# auxiliary functions
def find_pathway(graph, source, target):
    try:
        path = nx.shortest_path(graph, source=source, target=target, weight="distance")
        length = nx.shortest_path_length(graph, source=source, target=target, weight="distance")
    except nx.NetworkXNoPath:
        path, length = [], float('inf')
    return (source, target, path, length)

def pathway_strength(path, df_edges):
    """
    计算路径 path 的平均边权重
    path: list of nodes, e.g., [10,12,20,35]
    df_edges: DataFrame with columns ['Resid1','Resid2','weight']
    """
    weights = []
    for i in range(len(path)-1):
        u, v = path[i], path[i+1]
        
        # 因为边可能无向，所以查 u-v 或 v-u
        mask = ((df_edges['Resid1']==u) & (df_edges['Resid2']==v)) | \
               ((df_edges['Resid1']==v) & (df_edges['Resid2']==u))
        edge_weight = df_edges.loc[mask, 'weight'].values
        if len(edge_weight) > 0:
            weights.append(edge_weight[0])
        else:
            weights.append(0)  # 没有边时可以置0或nan
        
    if weights:
        return sum(weights)/len(weights)
    else:
        return 0

def ensure_list(x):
    if isinstance(x, list):
        return x
    elif isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
        return list(x)
    else:
        return [x]
    
def get_pairs(list1, list2):
    pairs = list(itertools.product(list1, list2))
    return pairs

# Application function
def app_pathway_between_regions(graph, df_edges, region1, region2):
    region1 = ensure_list(region1)
    region2 = ensure_list(region2)
    list_pairs = get_pairs(region1, region2)
    records = []
    for source, target in list_pairs:
        record = find_pathway(graph, source, target)
        records.append(record)
        
    df_paths = pd.DataFrame(records, columns=["Resid1", "Resid2", "path", "length"])
    df_paths.loc[:, 'strength'] = df_paths['path'].apply(lambda node: pathway_strength(node, df_edges))
    return df_paths


def app_pathway_to_distant(graph, df_edges, df_distance, region1, dist_threshold):
    region1 = ensure_list(region1)
    list_pairs = []
    df = df_distance.copy()
    df['pair'] = list(zip(df['Resid1'], df['Resid2']))
    for resid in region1:
        df_filtered = df[(df['Resid1'] == resid) | (df['Resid2'] == resid)]
        df_filtered = df_filtered[df_filtered['distance'] > float(dist_threshold)]
        list_pairs.extend(df_filtered["pair"].tolist())
    
    # list_pairs = list(set(list_pairs))
    records = []
    for source, target in list_pairs:
        record = find_pathway(graph, source, target)
        records.append(record)
        
    df_paths = pd.DataFrame(records, columns=["Resid1", "Resid2", "path", "length"])
    df_paths.loc[:, 'strength'] = df_paths['path'].apply(lambda node: pathway_strength(node, df_edges))
    return df_paths
    
    
def app_pathway_between_DccmClusters(graph, df_edges, cluster_DCCM, cluster1, cluster2):
    cluster1 = int(cluster1)
    cluster2 = int(cluster2)
    region1 = cluster_DCCM[cluster1]
    region2 = cluster_DCCM[cluster2]
    df_paths = app_pathway_between_regions(graph, df_edges, region1, region2)
    return df_paths

def app_cluster_DCCM(df_dccm, num_cluster, threshold=0.5):
    # 使用阈值挑选强协同节点
    dccm = df_dccm.values.copy()
    threshold = float(threshold)  # 可调整
    adj_matrix = dccm.copy()
    adj_matrix[np.abs(adj_matrix) < threshold] = 0  # 低相关置 0

    clustering = SpectralClustering(
        n_clusters=num_cluster, affinity='precomputed',
        assign_labels='kmeans', random_state=42
    ).fit(adj_matrix)
    labels = clustering.labels_
    
    dict_DCCM_clusters = {}
    # ========== 8. 保存每个 cluster 的残基编号 ==========
    for k in range(num_cluster):
        members = np.where(labels == k)[0].tolist()
        dict_DCCM_clusters[k+1] = [member + 1 for member in members]

    # ========== 9. 计算 cluster–cluster 平均相关性 ==========
    cluster_corr = np.zeros((num_cluster, num_cluster))
    for i in range(num_cluster):
        for j in range(num_cluster):
            members_i = np.where(labels == i)[0]
            members_j = np.where(labels == j)[0]
            sub_corr = dccm[np.ix_(members_i, members_j)]
            cluster_corr[i, j] = np.mean(sub_corr)

    # np.save("cluster_corr_matrix.npy", cluster_corr)
    df_cluster_corr = pd.DataFrame(cluster_corr)
    return dict_DCCM_clusters, df_cluster_corr

def app_dominant_motion(universe, df_pca, out_dir, target="PC1", rank=0.05):
    os.makedirs(out_dir, exist_ok=True)
    u = universe.copy()
    align.AlignTraj(u, u, select="protein and name CA", in_memory=True).run()
    if isinstance(target, list):
        target = target[0]
    df_sorted = df_pca[[target, "frame"]]
    df_sorted = df_sorted.sort_values(target).reset_index(drop=True)
    # 计算对应索引位置（5%、50%、95%）
    n = len(df_sorted)
    if rank > 0.5 and rank < 1:
        rank = rank - 0.5
    assert (rank > 0 and rank < 0.5)
    rank1 = float(rank) / 100
    rank2 = 1.0 - rank1
    positions = [int(n * rank1), int(n * 0.5), int(n * rank2)]
    # 防止越界（特别是样本量较小）
    positions = [min(max(p, 0), n - 1) for p in positions]
    # 取出对应行
    selected_rows = df_sorted.iloc[positions]
    frame_list = selected_rows["frame"].tolist()

    for i, frame in enumerate(frame_list):
        out_pdb = f"{out_dir}/{target}_{i}.pdb"
        with mda.Writer(out_pdb, multiframe=False) as W:
            u.trajectory[frame]
            W.write(u)
            
    return selected_rows
