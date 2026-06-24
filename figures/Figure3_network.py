#!/usr/bin/env python3
"""
Figure 3 — COP feature-based molecular network.

Renders the GNPS/MZmine IIMN network (connected clusters only) using the layout the
user built in Cytoscape (results.cys), with the Figure-1 group palette.
  - node = ion-identity-collapsed feature, labeled with its MZmine feature ID;
    pie = relative group contribution (heated 1000 / heated 100 / unheated 10 ug/mL);
    size = sum(precursor intensity); black outline = spectral-library match.
  - edge width = modified-cosine similarity.
  - clusters are numbered (Cluster 1..N) for reference in the text.
  - a PEG contaminant cluster is removed; singletons are not shown.

Self-contained. Run from anywhere (paths resolved relative to this file).
"""
import os, re, csv, zipfile
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe

# ---------------------------------------------------------------- paths
HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
GP   = os.path.join(BASE, "fbmn")
GRAPHML = os.path.join(GP, "iimn_collapsed.graphml")
CYS     = os.path.join(GP, "results.cys")
ANNOT   = os.path.join(BASE, "data", "output_annotations.csv")
OUT_PNG = os.path.join(HERE, "Fig3_network.png")
OUT_PDF = os.path.join(HERE, "Fig3_network.pdf")

# ---------------------------------------------------------------- palette (Figure 1)
C1000, C100, CUN = "#2CA58D", "#E8A33D", "#7E8B99"
GROUP_COLORS = [C1000, C100, CUN]
NODE_EDGE, NODE_EDGE_ANN = "#888888", "#111111"
EMPTY_FILL = "#D9D9D9"

# ID -> identity is reported in the figure caption (not drawn on the plot):
#   336 = 7α-hydroxy-4-cholesten-3-one, 488 = Cholest-4,6-dien-3-one,
#   522 = cholesterol, 1 = tetramethylpiperidine derivative
# Connected clusters seeded by these feature IDs are removed as contaminants:
REMOVE = {"PEG": "75", "polymer": "587"}   # PEG (tetraoxa-alcohol) + high-mass polymer series

# ---------------------------------------------------------------- parse graphml
G = "{http://graphml.graphdrawing.org/xmlns}"
root = ET.parse(GRAPHML).getroot()
keys = {k.get("id"): k.get("attr.name") for k in root.iter(G + "key")}
graph = root.find(G + "graph")
def node_data(n): return {keys[d.get("key")]: d.text for d in n.findall(G + "data")}

nodes = {}
for n in graph.findall(G + "node"):
    nodes[n.get("id")] = node_data(n)

edges = []                               # real edges only (drop GNPS self-loops)
for e in graph.findall(G + "edge"):
    s, t = e.get("source"), e.get("target")
    if s == t:
        continue
    cos = next((float(d.text) for d in e.findall(G + "data")
                if keys[d.get("key")] == "cosine_score"), 0.7)
    edges.append((s, t, cos))

# ---------------------------------------------------------------- layout from results.cys
zf = zipfile.ZipFile(CYS)
vname = [n for n in zf.namelist() if "/views/" in n and n.endswith(".xgmml")][0]
vx = zf.read(vname).decode("utf-8")
pos = {}
for m in re.finditer(
        r'<node id="\d+" label="([^"]+)"[^>]*>\s*<graphics z="[^"]*" x="([^"]+)" y="([^"]+)"', vx):
    pos[m.group(1)] = (float(m.group(2)), -float(m.group(3)))

# ---------------------------------------------------------------- annotations
ann = {}
with open(ANNOT) as f:
    for r in csv.DictReader(f):
        cid, sc = r["id"], float(r["score"] or 0)
        if cid not in ann or sc > ann[cid][1]:
            ann[cid] = (r["compound_name"].strip('"'), sc)

# ---------------------------------------------------------------- select + group
in_edge = set()
for s, t, _ in edges:
    in_edge.add(s); in_edge.add(t)
keep = [nid for nid, d in nodes.items() if nid in in_edge]
keep_set = set(keep)

parent = {nid: nid for nid in keep}
def find(a):
    while parent[a] != a:
        parent[a] = parent[parent[a]]; a = parent[a]
    return a
def union(a, b): parent[find(a)] = find(b)
for s, t, _ in edges:
    if s in keep_set and t in keep_set:
        union(s, t)
groups = {}
for nid in keep:
    groups.setdefault(find(nid), []).append(nid)

# drop contaminant clusters (whole connected components seeded by REMOVE)
ci2id = {nodes[n].get("cluster index"): n for n in keep}
removed = {}
for name, seed in REMOVE.items():
    mem = []
    if seed in ci2id:
        mem = groups.pop(find(ci2id[seed]), [])
        for nid in mem:
            keep_set.discard(nid)
    removed[name] = mem
keep = [n for n in keep if n in keep_set]

# ---------------------------------------------------------------- sizes
def fnum(d, k, default=0.0):
    try: return float(d.get(k))
    except (TypeError, ValueError): return default
XP, DIA = [0, 4, 5, 6, 7, 8], [10, 15, 22, 31, 42, 54]
def radius(intensity): return float(np.interp(np.log10(max(intensity, 1.0)), XP, DIA)) / 2.0

# ---------------------------------------------------------------- per-cluster layout
def relax(coords, rad, iters=160, pad=5.0):
    ks = list(coords); P = {k: np.array(coords[k], float) for k in ks}
    for _ in range(iters):
        for i in range(len(ks)):
            for j in range(i + 1, len(ks)):
                a, b = ks[i], ks[j]
                d = P[b] - P[a]; dist = float(np.hypot(*d)); need = rad[a] + rad[b] + pad
                if dist < 1e-6:
                    P[b] += np.array([need * 0.5, 0.0]); continue
                if dist < need:
                    u = d / dist; push = (need - dist) * 0.5
                    P[a] -= u * push; P[b] += u * push
    return {k: (float(P[k][0]), float(P[k][1])) for k in ks}

TARGET = 62.0
group_info = {}
for gkey, members in groups.items():
    xs = [pos[nodes[n].get("cluster index")][0] for n in members]
    ys = [pos[nodes[n].get("cluster index")][1] for n in members]
    rad = {n: radius(fnum(nodes[n], "sum(precursor intensity)")) for n in members}
    cx, cy = np.mean(xs), np.mean(ys)
    dd = max((np.hypot(xs[i] - cx, ys[i] - cy)) for i in range(len(members)))
    s = 1.0 if dd < 1e-6 else min(max(TARGET / dd, 0.05), 3.0)
    local = {n: ((xs[i] - cx) * s, (ys[i] - cy) * s) for i, n in enumerate(members)}
    local = relax(local, rad)                       # remove node overlaps
    cx2 = np.mean([p[0] for p in local.values()]); cy2 = np.mean([p[1] for p in local.values()])
    local = {n: (p[0] - cx2, p[1] - cy2) for n, p in local.items()}
    extent = max(np.hypot(*local[n]) + rad[n] for n in members)
    group_info[gkey] = dict(members=members, local=local, rad=rad, extent=max(extent, 18.0))

order = sorted(group_info, key=lambda g: (-len(group_info[g]["members"]), g))
cluster_no = {g: i + 1 for i, g in enumerate(order)}     # Cluster 1..N (largest first)

max_ext = max(gi["extent"] for gi in group_info.values())
cellw = 2 * max_ext + 40
cellh = 2 * max_ext + 46                                  # room for "Cluster N"
ncols = 3 if len(order) <= 6 else 4
node_xy, group_center = {}, {}
for idx, gkey in enumerate(order):
    r, c = divmod(idx, ncols)
    Cx, Cy = c * cellw, -r * cellh
    group_center[gkey] = (Cx, Cy)
    for nid, (lx, ly) in group_info[gkey]["local"].items():
        node_xy[nid] = (Cx + lx, Cy + ly)

# ---------------------------------------------------------------- draw
fig, ax = plt.subplots(figsize=(12, 8.5))
ax.set_aspect("equal"); ax.axis("off")

for s, t, cos in edges:                              # edges
    if s in keep_set and t in keep_set:
        x1, y1 = node_xy[s]; x2, y2 = node_xy[t]
        ax.plot([x1, x2], [y1, y2], "-", color="#CFCFCF", alpha=0.6,
                lw=float(np.interp(cos, [0.6, 0.8, 1.0], [0.5, 1.2, 2.4])),
                zorder=1, solid_capstyle="round")

for nid in keep:                                     # nodes
    d = nodes[nid]; ci = d.get("cluster index")
    x, y = node_xy[nid]; r = group_info[find(nid)]["rad"][nid] if find(nid) in group_info else radius(0)
    r = radius(fnum(d, "sum(precursor intensity)"))
    vals = [fnum(d, "GNPSGROUP:heated_1000ug"), fnum(d, "GNPSGROUP:heated_100ug"),
            fnum(d, "GNPSGROUP:unheated_10ug")]
    tot = sum(vals); is_ann = ci in ann
    ew, ec = (1.6, NODE_EDGE_ANN) if is_ann else (0.6, NODE_EDGE)
    if tot <= 0:
        ax.add_patch(Circle((x, y), r, facecolor=EMPTY_FILL, edgecolor=ec, lw=ew, zorder=3))
    else:
        a0 = 90.0
        for v, col in zip(vals, GROUP_COLORS):
            if v <= 0: continue
            a1 = a0 + 360.0 * v / tot
            ax.add_patch(Wedge((x, y), r, a0, a1, facecolor=col, edgecolor="none", zorder=3))
            a0 = a1
        ax.add_patch(Circle((x, y), r, facecolor="none", edgecolor=ec, lw=ew, zorder=4))
    ax.text(x, y, ci, ha="center", va="center", fontsize=5.6, color="#111111", zorder=7,
            path_effects=[pe.withStroke(linewidth=1.4, foreground="white")])

# "Cluster N" titles, aligned to a fixed height above each grid row
title_off = max_ext + 16
for gkey in order:
    cx, cy = group_center[gkey]
    ax.text(cx, cy + title_off, f"Cluster {cluster_no[gkey]}", ha="center", va="bottom",
            fontsize=9.5, fontweight="bold", color="#333333")

# limits (tight)
allx = [p[0] for p in node_xy.values()]; ally = [p[1] for p in node_xy.values()]
mx = max_ext + 22
ax.set_xlim(min(allx) - mx, max(allx) + mx)
ax.set_ylim(min(ally) - mx, max(ally) + title_off + 14)

# ---------------------------------------------------------------- legends (vertical),
# placed in the empty bottom-right grid cell to cut whitespace
ncells = ((len(order) - 1) // ncols + 1) * ncols
empty_idx = len(order) if len(order) < ncells else None
if empty_idx is not None:
    er, ec = divmod(empty_idx, ncols)
    ex, ey = ec * cellw, -er * cellh
else:                                          # fallback: just right of the grid
    ex, ey = (ncols) * cellw, 0.0

grp = [Line2D([0], [0], marker="o", linestyle="none", markersize=11,
              markerfacecolor=c, markeredgecolor="#555", label=l)
       for c, l in [(C1000, "Heated 1000 µg/mL"), (C100, "Heated 100 µg/mL"),
                    (CUN, "Unheated 10 µg/mL")]]
grp.append(Line2D([0], [0], marker="o", linestyle="none", markersize=11,
                  markerfacecolor="none", markeredgecolor=NODE_EDGE_ANN, markeredgewidth=1.8,
                  label="Library match"))
leg1 = ax.legend(handles=grp, loc="upper left",
                 bbox_to_anchor=(ex - max_ext, ey + max_ext), bbox_transform=ax.transData,
                 ncol=1, frameon=False, fontsize=9, handletextpad=0.5, labelspacing=0.9)
ax.add_artist(leg1)
size_handles = [Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="#cccccc",
                       markeredgecolor="#555", markersize=radius(v) * 0.55, label=f"{v:,.0f}")
                for v in (1e4, 1e5, 1e6)]
ax.legend(handles=size_handles, loc="upper left",
          bbox_to_anchor=(ex - max_ext, ey - max_ext * 0.18), bbox_transform=ax.transData,
          ncol=1, frameon=False, fontsize=8, title="Sum precursor intensity",
          title_fontsize=9, labelspacing=1.5, handletextpad=0.8)

fig.subplots_adjust(left=0.02, right=0.98, top=0.97, bottom=0.02)
fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
fig.savefig(OUT_PDF, bbox_inches="tight")

# ---------------------------------------------------------------- report
n_edges = sum(1 for s, t, _ in edges if s in keep_set and t in keep_set)
print(f"kept nodes      : {len(keep)} in {len(order)} clusters")
print(f"kept edges      : {n_edges}")
for name, mem in removed.items():
    cis = sorted((nodes[n].get('cluster index') for n in mem), key=lambda z: int(z))
    print(f"removed {name:>8} : {len(mem)} nodes ({cis})")
print("cluster numbering (Cluster N : #nodes : annotated CIs):")
for g in order:
    mem = group_info[g]["members"]
    anns = [nodes[n].get('cluster index') for n in mem
            if nodes[n].get('cluster index') in ann or nodes[n].get('cluster index') == '522']
    print(f"  Cluster {cluster_no[g]:>2} : {len(mem):>2} nodes : {sorted(anns, key=lambda z:int(z))}")
print(f"saved           : {OUT_PNG}")
