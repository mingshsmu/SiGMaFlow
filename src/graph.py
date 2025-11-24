import networkx as nx
import pandas as pd
import numpy as np
import MDAnalysis as mda
from itertools import combinations


def build_graph(
    universe: mda.Universe,
    df_ref_pdb: pd.DataFrame,
    df_weight: pd.DataFrame = None,
    eps = 1e-9
):
    """
    构建残基图，边权为 freq 或 freq*NMI，节点带 Reslabel 和 Category。
    
    Args:
        df_ref_pdb (pd.DataFrame): 残基信息表，必须包含 Resid, Resname, Category, Reslabel
    Returns:
        graph (nx.Graph): NetworkX 图，节点属性：resname, category, reslabel；边权 weight
        df_edges (pd.DataFrame): 边列表
    """
    u = universe.copy()
    # 选残基
    residues = [res for res in u.residues if res.resid in df_ref_pdb['Resid'].values]
    # if end is not None:
    #     residues = [res for res in residues if res.resid <= end]
    
    n_res = len(residues)
    resid_index = {res.resid: i for i, res in enumerate(residues)}
    
    # 创建图，添加节点属性
    graph = nx.Graph()
    for res in residues:
        ref_row = df_ref_pdb[df_ref_pdb['Resid'] == res.resid].iloc[0]
        graph.add_node(
            res.resid,
            resname = ref_row['Resname'],
            category = ref_row['Category'],
            reslabel = ref_row['Reslabel'],
            resid = ref_row['Resid']
        )
    
    # 构建边
    edges_data = []
    for i, j in combinations(range(n_res), 2):
        if df_weight is not None:
            try:
                weight = df_weight.at[residues[i].resid, residues[j].resid]
                if pd.isna(weight):
                    weight = 0.0
            except KeyError:
                weight = 0.0

        distance = 1.0 / (weight + eps)
        graph.add_edge(residues[i].resid, residues[j].resid, weight=weight, distance=distance)
        edges_data.append({"Resid1": residues[i].resid, "Resid2": residues[j].resid, "weight": weight, "distance": distance})
    
    df_edges = pd.DataFrame(edges_data)
    
    # node importance
    # 各种中心性指标
    degree_centrality = nx.degree_centrality(graph)
    betweenness_centrality = nx.betweenness_centrality(graph, weight="weight", normalized=True)
    closeness_centrality = nx.closeness_centrality(graph, distance="distance")
    try:
        eigen_centrality = nx.eigenvector_centrality(graph, weight="weight", max_iter=500)
    except nx.PowerIterationFailedConvergence:
        eigen_centrality = {n: 0 for n in graph.nodes()}
        print("! Eigenvector centrality did not converge; set to 0.")

    # Strength（加权度）
    strength = {n: sum(d.get('weight', 0.0) for _, _, d in graph.edges(n, data=True))
                for n in graph.nodes()}

    for n in graph.nodes():
        graph.nodes[n]["Degree"] = degree_centrality[n]
        graph.nodes[n]["Strength"] = strength[n]
        graph.nodes[n]["Betweenness"] = betweenness_centrality[n]
        graph.nodes[n]["Closeness"] = closeness_centrality[n]
        graph.nodes[n]["Eigenvector"] = eigen_centrality[n]

    # 整合结果
    df_importance = pd.DataFrame({
        'Resid': list(graph.nodes()),
        'Degree': [degree_centrality[n] for n in graph.nodes()],
        'Strength': [strength[n] for n in graph.nodes()],
        'Betweenness': [betweenness_centrality[n] for n in graph.nodes()],
        'Closeness': [closeness_centrality[n] for n in graph.nodes()],
        'Eigenvector': [eigen_centrality[n] for n in graph.nodes()]
    })

    # 排序
    df_importance = pd.merge(df_importance, df_ref_pdb[['Resid', 'Resname', 'Category', 'Reslabel']], on='Resid', how='left')
    df_importance.sort_values("Strength", ascending=False, inplace=True)
    df_importance.reset_index(drop=True, inplace=True)
    
    return graph, df_edges, df_importance

