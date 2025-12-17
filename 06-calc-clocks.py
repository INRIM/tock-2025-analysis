import decimal
import os
import time

import matplotlib.pyplot as plt
import numpy as np
import tintervals as ti
import tintervals.rocitlinks as rl

plt.close("all")
plt.ioff()

t0 = time.time()

dir = "./Data"
outdir = "./Analysis output"
figdir = os.path.join(outdir, "Figures")

os.makedirs(outdir, exist_ok=True)
os.makedirs(figdir, exist_ok=True)


# CCTF 2021 recommendations without rounding!
# used as a self-consistent reference

v0_fp = {
    "IT-Yb1": decimal.Decimal("518295836590863.630494915"),
    # "PTB-In+": decimal.Decimal("1267402452901041.28250601"),
    "PTB-In1": decimal.Decimal("1267402452901041.28250601"),
    "PTB-Al+": decimal.Decimal("1121015393207859.15926219"),
    "PTB-Yb1E3": decimal.Decimal("642121496772645.118522185"),
    "PTB-Sr3": decimal.Decimal("429228004229872.992467107"),
    "PTB-Sr4": decimal.Decimal("429228004229872.992467107"),
    # "SYRTE-Hg": decimal.Decimal("1128575290808154.31910196"),
    "OBSPARIS-Sr2": decimal.Decimal("429228004229872.992467107"),
    "OBSPARIS-SrB": decimal.Decimal("429228004229872.992467107"),
    # "PTB-Yb1E3E2": decimal.Decimal("688358979309308.239120953"),
    "NPL-E3Yb+3": decimal.Decimal("642121496772645.118522185"),
    "NPL-Sr1": decimal.Decimal("429228004229872.992467107"),
    "INRIM_HM": decimal.Decimal("1"),
    "PTB_HM": decimal.Decimal("1"),
}

long_names = {
    "IT-Yb1": "INRIM_ITYb1",
    "PTB-In1": "PTB_In_CombKnoten",
    "PTB-Al+": "PTB_Al_CombAl",
    "PTB-Yb1E3": "PTB_Yb_CombKnoten",
    "PTB-Sr3": "PTB_Sr3_CombKnoten",
    "PTB-Sr4": "INRIM_PTBSr4",
    "OBSPARIS-Sr2": "OBSPARIS_Sr2",
    "OBSPARIS-SrB": "OBSPARIS_SrB",
    "PTB-Yb1E3E2": "PTB_Yb1E2_CombYb",
    "NPL-E3Yb+3": "NPL_YbE3",
    "NPL-Sr1": "NPL_Sr1",
}


clocks = [
#    "PTB-In1",
#    "PTB-Al+",
    #"NPL-E3Yb+3",
#    "PTB-Yb1E3",
#    "IT-Yb1",
    "NPL-Sr1",
    "OBSPARIS-Sr2",
    "OBSPARIS-SrB",
    "PTB-Sr3",
#    "PTB-Sr4",
    ]


print("TOCK March 2025")

comparator_names = [f.name for f in os.scandir(dir) if f.is_dir()]


start = 60740
stop = 60900


links = {}
connections = []

print("Loading files")
for name in comparator_names:
    # tstart, tstop = startstop[name]
    tstart, tstop = start, stop
    try:
        links[name] = rl.load_link_from_dir(
            os.path.join(dir, name), start=ti.epoch_from_mjd(tstart), stop=ti.epoch_from_mjd(tstop)
        )
    except IOError:
        print(f"cannot load {name}")
        continue
    print(f"loading {name}")
    connections += [name]

# NPL uncertainty is given in Hz, not in relative
# TODO: make code smarter
#links["NPL_T1-NPL_Sr1"].data[:, 3] /= float(links["NPL_T1-NPL_Sr1"].oscA.v0)
#links["NPL_T1-NPL_YbE3"].data[:, 3] /= float(links["NPL_T1-NPL_YbE3"].oscA.v0)

# # temporary fix for LPL_MLS-LPL_RLS timing offset
# links["LPL_MLS-LPL_RLS"].t += 3600.


# plot uptimes
fig, ax = plt.subplots(figsize=(10, 10))
ylab = []
N = len(links.values())


def sort_key(link):
    name = link.name
    # First, check if name contains any clock long name
    is_clock = any(long_names[sub] in name for sub in clocks)
    # Sort clock name first (False < True), then alphabetically
    return (not is_clock, name)


# fig.suptitle('Uptimes')
for i, link in enumerate(sorted(links.values(), key=sort_key, reverse=True)):
    vals = ti.array2intervals(link.t, tgap=100.0, tblock=100.0)
    vals = ti.mjd_from_epoch(vals)

    height = 0.8
    y_pos = i - height / 2
    for a, b in vals:
        ax.broken_barh([(a, b - a)], (y_pos, height), color="C{}".format(i))

    ylab += [link.name]


ax.set_xlim(start, stop)
ax.set_xlabel("MJD")
ax.set_ylim(-0.5, N - 0.5)
ax.set_yticks(np.arange(N), ylab)

ax.grid(axis="x", alpha=0.7)
fig.tight_layout()
fig.savefig(os.path.join(figdir, "network-uptimes.png"))


# graphing
def build_graph(connections):
    graph = {}
    for conn in connections:
        a, b = conn.split("-")
        graph.setdefault(a, []).append(b)
        graph.setdefault(b, []).append(a)
    return graph


def shortest_path(graph, start, end):
    from collections import deque

    if start not in graph or end not in graph:
        return None

    queue = deque([(start, [start])])
    visited = set()

    while queue:
        current, path = queue.popleft()
        if current == end:
            return [f"{path[i]}-{path[i+1]}" for i in range(len(path) - 1)]
        visited.add(current)
        for neighbor in graph.get(current, []):
            if neighbor not in visited and neighbor not in path:
                queue.append((neighbor, path + [neighbor]))
    return None


def invert_connection(conn):
    a, b = conn.split("-")
    return f"{b}-{a}"


graph = build_graph(connections)


chains = {}
print("Calculating chains")
for i, A in enumerate(clocks):
    for j, B in enumerate(clocks[i + 1 :]):
        short = shortest_path(graph, long_names[A], long_names[B])
        chains[A + "/" + B] = short
        print(A + "/" + B, short)

# add masers
chains["INRIM_HM/PTB_HM"] = shortest_path(graph, "INRIM_HM", "PTB_HM_CombKnoten")


# prepare all chains
ratios = {}
print("Building chains")
for name, comps in chains.items():
    print("\n", name)
    if comps:
        res = []
        for c in comps[::-1]:
            if c in (connections):
                print(c)
                res += [links[c]]
            elif invert_connection(c) in (connections):
                print(invert_connection(c))
                res += [-links[invert_connection(c)]]
            else:
                print("MISSING!:", c)

        ratios[name] = res


reslinks = {}

i = 0
print("Calculating clock ratios")
for shortname, llinks in list(ratios.items()):
    print("\n" + shortname)
    clock2, clock1 = shortname.split("/")

    if len(llinks) == 1:
        reslink = llinks[0]
    else:
        reslink, masks = rl.chain(*llinks)
    print(reslink.name)
    # for some reasons these prints are affected by later changes to the class

    # normalize to my choiche of r0
    reslink.oscA.v0 = v0_fp[clock1]
    reslink.oscB.v0 = v0_fp[clock2]
    reslink.normalize()
    print(reslink.r0)

    reslinks[shortname] = reslink

    # time varying uB
    if llinks[0].data.shape[1] == 4:
        uuB1 = llinks[0].data[masks[0], 3]
        # reslink.oscA.systematic_uncertainty = np.mean(uuB1)
    else:
        uuB1 = reslink.oscA.systematic_uncertainty * np.ones_like(reslink.t)

    # time varying uB
    if llinks[-1].data.shape[1] == 4:
        uuB2 = llinks[-1].data[masks[-1], 3]
        # reslink.oscB.systematic_uncertainty = np.mean(uuB2)
    else:
        uuB2 = reslink.oscB.systematic_uncertainty * np.ones_like(reslink.t)

    reslink.data = np.column_stack((reslink.data, uuB1, uuB2))

    rl.save_link_to_dir(outdir, reslink, extra_names=["usys_A", "usys_B"])


# save link name/shortname
shortnames = list(reslinks.keys())
shortnames.remove('INRIM_HM/PTB_HM')
longnames = [reslinks[x].name for x in shortnames]
names = np.column_stack((shortnames, longnames))
np.savetxt(os.path.join(outdir, "ratios.txt"), names, fmt="%s\t%s")


print(f"Saving clock data time = {time.time() - t0} s")
