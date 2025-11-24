# get_atom_coord_via_atomspec.py
from chimerax.core.commands.atomspec import AtomSpecArg
# atom_spec 示例可以是 "#1/A:25@CA" 或 ":25.A@CA" 等
atom_spec = "#1/X:66,72@CA"

# 解析 -> 返回 AtomSpecArg 对象列表（可能包含多个 selector）
parsed = AtomSpecArg.parse(atom_spec, session)
if not parsed:
    print("无法解析 atomspec:", atom_spec)
else:
    # 取第一个解析结果并 evaluate(session) -> 得到一个 AtomSpec 对象
    atom_spec_obj = parsed[0].evaluate(session)
    atoms = getattr(atom_spec_obj, "atoms", None)
    if not atoms:
        print("在该 atomspec 下没有匹配到原子。")
    else:
        atom = atoms[0]   # 第一个匹配的原子
        coords = [tuple(atom.coord) for atom in atoms]   # (x, y, z)
        residues = [f"{atom.residue.chain_id}:{atom.residue.number}" for atom in atoms]

        for res, c in zip(residues, coords):
            print(f"{res} -> ({c[0]:.3f}, {c[1]:.3f}, {c[2]:.3f})")
            

from chimerax.core.commands import run

# 定义一组圆锥的参数
cones = [
    {"fromPoint": (0, 0, 0), "toPoint": (0, 0, 10), "topRadius": 0.1, "radius": 3, "color": "red"},
    {"fromPoint": (0, 0, 0), "toPoint": (0, 10, 0), "topRadius": 0.1, "radius": 3, "color": "blue"},
    {"fromPoint": (0, 0, 0), "toPoint": (10, 0, 0), "topRadius": 0.1, "radius": 3, "color": "green"},
]

# 循环绘制
for c in cones:
    cmd = (
        f"shape cone fromPoint {c['fromPoint'][0]},{c['fromPoint'][1]},{c['fromPoint'][2]} "
        f"toPoint {c['toPoint'][0]},{c['toPoint'][1]},{c['toPoint'][2]} "
        f"radius {c['radius']} topRadius {c['topRadius']} color {c['color']}"
    )
    run(session, cmd)
