import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import networkx as nx
from scipy.stats import gaussian_kde
from matplotlib.patches import Patch

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


def plot_line(df, x, y, xlab, ylab, title, color="#fc8d62", save_path=None, fig_size=(8, 5)):
    plt.figure(figsize=fig_size)
    sns.lineplot(data=df, x=x, y=y, color=color, errorbar=None)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    if save_path is not None:
        save_path, _ = os.path.splitext(save_path)
        plt.savefig(f"{save_path}.png", dpi=600)
        plt.savefig(f"{save_path}.pdf")
    
    plt.show()


def plot_heatmap(df, xlab, ylab, vmin, vmax, title, colorbar_label, axis_break=None,cmap='bwr', origin='lower', save_path=None, fig_size=(6, 5)):
    if axis_break is None:
        list_break = [1, 2, 5, 10, 20, 50, 100]
        for b in list_break:
            ticks_count = df.shape[0] / b + 1
            if ticks_count > 4 and ticks_count < 10:
                axis_break = b
                break
    
    range_ticks = range(0, df.shape[0], axis_break)
    label_ticks = [i+1 for i in range(0, df.shape[0], axis_break)]
    plt.figure(figsize=fig_size)
    plt.imshow(df.values, cmap=cmap, vmin=vmin, vmax=vmax, origin=origin)
    plt.colorbar(label=colorbar_label)
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.xticks(range_ticks, label_ticks)
    plt.yticks(range_ticks, label_ticks)
    plt.tight_layout()
    if save_path is not None:
        save_path, _ = os.path.splitext(save_path)
        plt.savefig(f"{save_path}.png", dpi=600)
        plt.savefig(f"{save_path}.pdf")
    plt.show()


def plot_scatter(df, x, y, xlab, ylab, title, cmap='viridis', 
                 color_strategy=['time_series', 'density'], save_path=None, fig_size=(6, 5)):
    if isinstance(color_strategy, list):
        color_strategy = color_strategy[0]
    
    x_np = df[x].values
    y_np = df[y].values
    if color_strategy == "time_series":
        colorbar_label = "Frame index"
        z = df.index
        
    elif color_strategy == "density":
        colorbar_label = "Density"
        xy = np.vstack([x_np, y_np])
        z = gaussian_kde(xy)(xy)
        
    plt.figure(figsize=fig_size)
    plt.scatter(x_np, y_np, s=10, c=z, cmap=cmap)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    plt.colorbar(label=colorbar_label)
    if save_path is not None:
        save_path, _ = os.path.splitext(save_path)
        plt.savefig(f"{save_path}.png", dpi=600)
        plt.savefig(f"{save_path}.pdf")
    plt.show()
    
    
def plot_graph(graph, fig_size=(10, 8), category_colors=None, save_path=None, font_size=6,
               node_label=["reslabel", "resid", "resname"]):
    pos = nx.spring_layout(graph, seed=42, k=0.1)
    # pos = nx.kamada_kawai_layout(graph_residue)
    # pos = nx.spectral_layout(graph_residue)
    # pos = nx.circular_layout(graph_residue)
    # pos = nx.shell_layout(graph_residue)
    # 提取节点属性
    if isinstance(node_label, list):
        node_label = node_label[0]
    
    with_label = node_label is not None
    
    if with_label:
        if node_label in ["reslabel", "resid", "resname"]:
            labels = nx.get_node_attributes(graph, node_label)
        else:
            raise ValueError(f"Unsupported node_label: {node_label}")
    
    categories = nx.get_node_attributes(graph, "category")

    # 定义颜色映射（根据类别）
    if category_colors is None:
        category_colors = {
            "Acidic": "#ed028c",
            "Basic": "#40c7f4",
            "Polar": "#fad5e5",
            "Hydrophobic": "#b9e5fa",
            "Aromatic": "#f5c45e",
            "Special": "#a5d399",
            "Neutral": "#d9d9d9",
        }

    node_colors = [category_colors.get(categories.get(n, "Neutral"), "gray")
                for n in graph.nodes()]

    # 可根据边权调整宽度
    widths = [d["weight"] * 5 for _, _, d in graph.edges(data=True)]

    # 节点重要性 Size
    # strength_dict = dict(zip(df_importance["Resid"], df_importance["Strength"]))
    # nx.set_node_attributes(graph, strength_dict, "Strength")
    # strengths = np.array([strength_dict.get(n, 0) for n in graph_residue.nodes()])
    # # 归一化到 [min_size, max_size]
    min_size, max_size = 50, 1000  # 你可以调节范围
    strength_dict = nx.get_node_attributes(graph, "Strength")
    strengths = np.array([strength_dict.get(node, 0) for node in graph.nodes()])

    # 排名（越大排名越靠后）
    ranks = np.argsort(np.argsort(strengths))

    # 将排名线性映射到 [min_size, max_size]
    if len(ranks) > 1:
        node_sizes = min_size + (ranks / (len(ranks) - 1)) * (max_size - min_size)
    else:
        node_sizes = np.array([min_size])


    # 绘图
    plt.figure(figsize=fig_size)
    if with_label:
        nx.draw(
            graph,
            pos,
            with_labels=True,
            labels=labels,              # 自定义 label
            node_color=node_colors,     # 每个节点不同颜色
            edge_color='gray',
            width=widths,
            node_size=node_sizes,
            font_size=font_size,
            font_weight='bold')
    else:
        nx.draw(
            graph,
            pos,
            with_labels=False,
            node_color=node_colors,     # 每个节点不同颜色
            edge_color='gray',
            width=widths,
            node_size=node_sizes)


    legend_elements = [Patch(facecolor=color, label=cat)
                    for cat, color in category_colors.items()]
    plt.legend(handles=legend_elements, title="Residue Type")
    if save_path is not None:
        save_path, _ = os.path.splitext(save_path)
        plt.savefig(f"{save_path}.png", dpi=600)
        plt.savefig(f"{save_path}.pdf")
    plt.show()
    
    