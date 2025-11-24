import platform
import os

# to avoid Kmeans memory leak on windows
if platform.system() == "Windows":
    if "OMP_NUM_THREADS" not in os.environ: 
        os.environ["OMP_NUM_THREADS"] = "1"

from Bio.PDB import PDBParser
import pandas as pd
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

from tqdm import tqdm
from sklearn.cluster import SpectralClustering
from sklearn.metrics import silhouette_score, normalized_mutual_info_score

from multiprocessing import Pool, cpu_count

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

import MDAnalysis as mda
from MDAnalysis.analysis import pca, rms, align, diffusionmap
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.analysis.dihedrals import Dihedral

from collections.abc import Iterable
import itertools
from itertools import combinations

import re
import json
import logging

from .prepare import *
from .compute import *
from .graph import *
from .application import *
from .plot import *
from .visualize import *



logging.basicConfig(level=logging.INFO)

# Main class
class SiGMaFlow():
    def __init__(self, topology, trajectory, reference):
        self.topology = topology
        self.trajectory = trajectory
        self.reference = reference
    
    # ==========================
    # === Module: 1. prepare ===
    # ==========================
    def prepare(self):
        logging.info(f"[~] Loading trajectory from {self.trajectory} and processing reference structure from {self.reference}.")
        self.df_ref_pdb = parse_pdb(self.reference)
        self.df_AApair_dist = calculate_AA_dist(self.df_ref_pdb)
        self.u = mda.Universe(self.topology, self.trajectory)
    
    # ==========================
    # === Module: 2. compute ===
    # ==========================
    def computePhysicalQuantities(self, physical_quantities = ["Dihedrals", "RMSF", "PCA", "DCCM"]):
        if isinstance(physical_quantities, str):
            physical_quantities = [physical_quantities]
        for physical_quantity in physical_quantities:
            match physical_quantity:
                case "Dihedrals":
                    logging.info("[~] Computing Dihedrals... Result will store into self.df_dihedrals")
                    self.df_dihedrals = compute_dihedrals(universe=self.u)
                
                case "RMSF":
                    logging.info("[~] Computing RMSF... Result will store into self.df_rmsf")
                    self.df_rmsf = compute_RMSF(universe=self.u)
                
                case "PCA":
                    logging.info("[~] Conducting PCA... Result will store into self.df_pca")
                    self.df_pca = compute_PCA(universe=self.u)
                
                case "DCCM":
                    logging.info("[~] Computing DCCM... Result will store into self.df_dccm")
                    self.df_dccm = compute_DCCM(universe=self.u)
                            
                            
    def computeContactFrequency(self, dist: float=5.0, heavy_atoms: list=["C","N","O","S"]):
        logging.info("[~] Computing Residue contact frequency... Result will store into self.df_contact")
        self.df_contact = compute_contact_frequency(universe=self.u, df_ref=self.df_ref_pdb, dist=dist, heavy_atoms=heavy_atoms)
        
    def computeNMI(self, n_bins=72, average_method="geometric"):
        logging.info("[~] Computing Dihedrals NMI... Result will store into self.df_nmi_phi, self.df_nmi_psi, self.df_nmi_chi1")
        logging.info("[~] Phi...")
        self.df_nmi_phi = compute_NMI(df=self.df_dihedrals, obj="Phi", n_bins=n_bins, average_method=average_method)
        logging.info("[~] Psi...")
        self.df_nmi_psi = compute_NMI(df=self.df_dihedrals, obj="Psi", n_bins=n_bins, average_method=average_method)
        logging.info("[~] Chi1...")
        self.df_nmi_chi1 = compute_NMI(df=self.df_dihedrals, obj="Chi1", n_bins=n_bins, average_method=average_method)
    
    def computeIntegratedNMI(self, weights=[0.1, 0.1, 0.8]):
        self.df_nmi = compute_weighted_average(dfs=[self.df_nmi_phi, self.df_nmi_psi, self.df_nmi_chi1],
                                               weights=weights)
    
    def computeEdgeWeight(self):
        self.df_weight = compute_edge_weight(df_NMI=self.df_nmi, df_contact=self.df_contact)
    
    # ========================
    # === Module: 3. graph ===
    # ========================
    def graphBuild(self):
        self.graph, self.df_edges, self.df_node_importance = build_graph(universe=self.u, df_ref_pdb=self.df_ref_pdb, 
                                                                         df_weight=self.df_weight, eps=1e-9)
    
    # ==============================
    # === Module: 4. application ===
    # ==============================
    def appClusterDCCM(self, num_cluster, threshold=0.5):
        self.dict_dccm_clusters, self.df_dccm_corr = app_cluster_DCCM(df_dccm=self.df_dccm, num_cluster=num_cluster,
                                                                      threshold=threshold)
        return self.dict_dccm_clusters
    
    def appPathwayBtwRegions(self, region1, region2):
        df_paths = app_pathway_between_regions(graph=self.graph, df_edges=self.df_edges, region1=region1, region2=region2)
        return df_paths
    
    def appPathwayToDistant(self, region1, dist_threshold):
        df_paths = app_pathway_to_distant(graph=self.graph, df_edges=self.df_edges, df_distance=self.df_AApair_dist,
                                          region1=region1, dist_threshold=dist_threshold)
        return df_paths
    
    def appPathwayBtwDccmClusters(self, cluster1, cluster2):
        df_paths = app_pathway_between_DccmClusters(graph=self.graph, df_edges=self.df_edges,
                                                    cluster_DCCM=self.dict_dccm_clusters,
                                                    cluster1=cluster1, cluster2=cluster2)
        return df_paths
    
    def appDominantMotion(self, out_dir="PCA_motion", target="PC1", rank=0.05):
        logging.info(f"[~] PCA Dominant Motion, will export PDBs into directory: {out_dir}")
        df_selected_rows = app_dominant_motion(universe=self.u, df_pca=self.df_pca, out_dir=out_dir, target=target, rank=rank)
        return df_selected_rows
    
    # =======================
    # === Module: 5. plot ===
    # =======================
    def plotLine(self, obj, x=None, y=None, xlab=None, ylab=None, title=None, color="#fc8d62", save_path=None, fig_size=(8, 5)):
        if isinstance(obj, str):
            if obj == "rmsf" or obj == "RMSF":
                plot_line(df=self.df_rmsf, x="resid", y="rmsf", xlab="Residue Index", ylab="RMSF (Å)", title=None, 
                        color=color, save_path=save_path, fig_size=fig_size)
            else:
                raise ValueError(f"Unsupported str 'obj': {obj}")
            
        else:
            plot_line(df=obj, x=x, y=y, xlab=xlab, ylab=ylab, title=title, 
                        color=color, save_path=save_path, fig_size=fig_size)

    def plotHeatmap(self, obj, xlab=None, ylab=None, vmin=None, vmax=None, title=None, axis_break=None, colorbar_label=None, cmap='bwr', origin='lower', save_path=None, fig_size=(6, 5)):
        if isinstance(obj, str):
            if obj == "DCCM" or obj == "dccm":
                plot_heatmap(df=self.df_dccm, xlab="Residue Index", ylab="Residue Index", vmin=-1, vmax=1, 
                             title="DCCM (CA atoms)", axis_break=axis_break, colorbar_label="Correlation coefficient", cmap=cmap,
                             origin="lower", save_path=save_path, fig_size=fig_size)
            elif obj == "DCCM_corr" or obj == "dccm_corr":
                plot_heatmap(df=self.df_dccm_corr, xlab="Cluster Index", ylab="Cluster Index", vmin=-1, vmax=1, 
                             title="DCCM Cluster Correlation Map", axis_break=1, colorbar_label="Correlation coefficient", cmap=cmap,
                             origin="lower", save_path=save_path, fig_size=fig_size)
            elif obj == "NMI_Phi" or obj == "nmi_phi":
                mat = self.df_nmi_phi.values
                vmin = np.nanmin(mat)
                vmax = np.nanmax(mat)
                mat_without_max = np.where(mat == vmax, -np.inf, mat)
                second_max = np.nanmax(mat_without_max)
                plot_heatmap(df=self.df_nmi_phi, xlab="Residue Index", ylab="Residue Index", vmin=vmin, vmax=second_max, 
                             title="NMI: Phi", axis_break=axis_break, colorbar_label="Normalized Mutual Information", cmap=cmap,
                             origin="lower", save_path=save_path, fig_size=fig_size)
            elif obj == "NMI_Psi" or obj == "nmi_psi":
                mat = self.df_nmi_psi.values
                vmin = np.nanmin(mat)
                vmax = np.nanmax(mat)
                mat_without_max = np.where(mat == vmax, -np.inf, mat)
                second_max = np.nanmax(mat_without_max)
                plot_heatmap(df=self.df_nmi_psi, xlab="Residue Index", ylab="Residue Index", vmin=vmin, vmax=second_max, 
                             title="NMI: Psi", axis_break=axis_break, colorbar_label="Normalized Mutual Information", cmap=cmap,
                             origin="lower", save_path=save_path, fig_size=fig_size)
            elif obj == "NMI_Chi1" or obj == "nmi_chi1":
                mat = self.df_nmi_chi1.values
                vmin = np.nanmin(mat)
                vmax = np.nanmax(mat)
                mat_without_max = np.where(mat == vmax, -np.inf, mat)
                second_max = np.nanmax(mat_without_max)
                plot_heatmap(df=self.df_nmi_chi1, xlab="Residue Index", ylab="Residue Index", vmin=vmin, vmax=second_max, 
                             title="NMI: Chi1", axis_break=axis_break, colorbar_label="Normalized Mutual Information", cmap=cmap,
                             origin="lower", save_path=save_path, fig_size=fig_size)
            elif obj == "NMI" or obj == "nmi":
                mat = self.df_nmi.values
                vmin = np.nanmin(mat)
                vmax = np.nanmax(mat)
                mat_without_max = np.where(mat == vmax, -np.inf, mat)
                second_max = np.nanmax(mat_without_max)
                plot_heatmap(df=self.df_nmi, xlab="Residue Index", ylab="Residue Index", vmin=vmin, vmax=second_max, 
                             title="Integrated NMI", axis_break=axis_break, colorbar_label="Normalized Mutual Information", cmap=cmap,
                             origin="lower", save_path=save_path, fig_size=fig_size)
            elif obj == "contact" or obj == "contact_frequency":
                plot_heatmap(df=self.df_contact, xlab="Residue Index", ylab="Residue Index", vmin=0, vmax=1, 
                             title="Contact", axis_break=axis_break, colorbar_label="Frequency", cmap=cmap,
                             origin="lower", save_path=save_path, fig_size=fig_size)
            elif obj == "weight":
                mat = self.df_weight.values
                vmin = np.nanmin(mat)
                vmax = np.nanmax(mat)
                plot_heatmap(df=self.df_weight, xlab="Residue Index", ylab="Residue Index", vmin=vmin, vmax=vmax, 
                             title="Weight: NMI * Contact", axis_break=axis_break, colorbar_label="Strength", cmap=cmap,
                             origin="lower", save_path=save_path, fig_size=fig_size)
            else:
                raise ValueError(f"Unsupported str 'obj': {obj}")
                
        else:
            plot_heatmap(df=obj, xlab=xlab, ylab=ylab, vmin=vmin, vmax=vmax, 
                             title=title, axis_break=axis_break, colorbar_label=colorbar_label, cmap=cmap,
                             origin=origin, save_path=save_path, fig_size=fig_size)
        
    def plotScatter(self, obj, x=None, y=None, xlab=None, ylab=None, title=None, cmap='viridis', color_strategy=None,
                    save_path=None, fig_size=(6, 5)):
        if isinstance(obj, str):
            if obj == "PCA" or obj == "pca":
                if color_strategy is None:
                    color_strategy = "time_series"
                if isinstance(color_strategy, list):
                    color_strategy = color_strategy[0]
                plot_scatter(df=self.df_pca, x="PC1", y="PC2", xlab="PC1", ylab="PC2", title="PCA projection (CA atoms)",
                             cmap=cmap, color_strategy=color_strategy, save_path=save_path, fig_size=fig_size)
            else:
                raise ValueError(f"Unsupported str 'obj': {obj}")
        else:
            plot_scatter(df=obj, x=x, y=y, xlab=xlab, ylab=ylab, title=title,
                             cmap=cmap, color_strategy=color_strategy, save_path=save_path, fig_size=fig_size)
        
    def plotGraph(self, df_paths_filtered=None,fig_size=(10, 8), category_colors=None, 
                  save_path=None, node_label=["reslabel", "resid", "resname"]):
        graph = self.graph.copy()
        if df_paths_filtered is not None:
            list_nodes = df_paths_filtered["path"].explode().unique().tolist()
            graph = graph.subgraph(nodes=list_nodes).copy()

        plot_graph(graph=graph, fig_size=fig_size, category_colors=category_colors, save_path=save_path, node_label=node_label)
        
        
    # ============================
    # === Module: 6. visualize ===
    # ============================
    def visRmsfPymol(self, outfile, cmap='coolwarm', midpoint=None, half_range=None):
        visualize_RMSF_pymol(df_rmsf=self.df_rmsf, outfile=outfile, cmap=cmap, midpoint=midpoint, half_range=half_range)
    
    def visRmsfChimeraX(self, outfile, cmap='coolwarm', midpoint=None, half_range=None):
        visualize_RMSF_chimeraX(df_rmsf=self.df_rmsf, outfile=outfile, cmap=cmap, midpoint=midpoint, half_range=half_range)
    
    def visRmsfColorbar(self, outfile=None, cmap='coolwarm', midpoint=None, half_range=None, fig_size=(0.5, 8), n_break=7):
        visualize_RMSF_bar(outfile=outfile, cmap=cmap, midpoint=midpoint, half_range=half_range, fig_size=fig_size, n_break=n_break)
    
    def visDccmClusterPymol(self, outfile, colors=None):
        visualize_DccmCluster_pymol(dict_dccm_clusters=self.dict_dccm_clusters, outfile=outfile, colors=colors)
        
    def visDccmClusterChimeraX(self, outfile, colors=None):
        visualize_DccmCluster_chimeraX(dict_dccm_clusters=self.dict_dccm_clusters, outfile=outfile, colors=colors)
        
    def visPathPreprocess(self, df_paths, outfile, cluster_by="path_id", colors=None, radius_range=(0.1, 1.0)):
        dict_coords = visualize_path_preprocess(df_paths=df_paths, df_edges=self.df_edges, df_ref_pdb=self.df_ref_pdb,
                                                outfile=outfile, cluster_by=cluster_by, colors=colors, radius_range=radius_range)
        return dict_coords
    
    def visPathPymol(self, dict_coords, outfile=None, colors=None):
        visualize_path_pymol(dict_coords=dict_coords, outfile=outfile, colors=colors)
        
    def visPathChimeraX(self, dict_coords, outfile=None, colors=None):
        visualize_path_chimeraX(dict_coords=dict_coords, outfile=outfile, colors=colors)
        
    def visMotionPymol(self, obj1, obj2, outfile, low_threshold=1, high_threshold=3,
                           cone_color="orange", cone_sep=1, cone_radius_range=(0.05, 0.5)):
        visualize_motion_pymol(obj1=obj1, obj2=obj2, df_rmsf=self.df_rmsf, outfile=outfile, low_threshold=low_threshold,
                               high_threshold=high_threshold, cone_color=cone_color, cone_sep=cone_sep, cone_radius_range=cone_radius_range)
        
    def visMotionChimeraX(self, obj1, obj2, outfile, low_threshold=1, high_threshold=3,
                           cone_color="orange", cone_sep=1, cone_radius_range=(0.05, 0.5)):
        visualize_motion_chimeraX(obj1=obj1, obj2=obj2, df_rmsf=self.df_rmsf, outfile=outfile, low_threshold=low_threshold,
                               high_threshold=high_threshold, cone_color=cone_color, cone_sep=cone_sep, cone_radius_range=cone_radius_range)    