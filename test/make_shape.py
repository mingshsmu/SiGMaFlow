# 文件名：make_cylinders.py

from chimerax.core.commands import run

# 示例数据：每个 cluster 的 cylinder 坐标、半径、颜色
dict_coords = {
    "cluster1": [
        {"coord1": [0, 0, 0], "coord2": [10, 0, 0], "radius": 0.3, "color": [1, 0, 0]},
        {"coord1": [0, 0, 0], "coord2": [0, 10, 0], "radius": 0.3, "color": [0, 1, 0]},
    ],
    "cluster2": [
        {"coord1": [0, 0, 0], "coord2": [0, 0, 10], "radius": 0.3, "color": [0, 0, 1]},
        {"coord1": [10, 0, 0], "coord2": [10, 10, 0], "radius": 0.2, "color": [1, 1, 0]},
    ]
}

# 保存 cluster -> list of modelIDs 的映射
cluster_models = {}

for cluster_name, rows in dict_coords.items():
    cluster_models[cluster_name] = []
    for i, row in enumerate(rows):
        model_name = f"{cluster_name}_{i}"
        p1 = ",".join(str(x) for x in row["coord1"])
        p2 = ",".join(str(x) for x in row["coord2"])
        r = row["radius"]
        # ChimeraX color 要求 0-1 范围
        color_str = ",".join(str(c*100) for c in row["color"])
        cmd = f"shape cylinder fromPoint {p1} toPoint {p2} radius {r} color {color_str} name {model_name}"
        run(session, cmd)
        # 获取刚创建模型的 ID（假设模型 ID 自增，可以用最后一个模型）
        # 注意：ChimeraX 1.10.1 中 run 返回 None，获取 model ID 需要额外操作
        # 简单方法：每创建一个模型，使用最后一个模型 ID
        modelID = session.models.list()[-1].id
        cluster_models[cluster_name].append(modelID)

# 定义函数，批量显示/隐藏 cluster
def show_cluster(cluster_name):
    if cluster_name in cluster_models:
        ids = " ".join(f"#{mid}" for mid in cluster_models[cluster_name])
        run(session, f"show {ids} models")

def hide_cluster(cluster_name):
    if cluster_name in cluster_models:
        ids = " ".join(f"#{mid}" for mid in cluster_models[cluster_name])
        run(session, f"hide {ids} models")

# 使用示例：
# show_cluster("cluster1")
# hide_cluster("cluster2")
