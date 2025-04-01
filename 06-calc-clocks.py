import decimal
import os
import time

import numpy as np
import tintervals as ti
import tintervals.rocitlinks as rl

t0 = time.time()

dir = "./Data"
outdir = "./Analysis output"

os.makedirs(outdir, exist_ok=True)

# CCTF 2021 recommendations without rounding!
# used as a self-consistent reference

v0_fp = {
    "IT-Yb1": decimal.Decimal("518295836590863.630494915"),
    # "PTB-In+": decimal.Decimal("1267402452901041.28250601"),
    "PTB-Sr3": decimal.Decimal("429228004229872.992467107"),
    "PTB-Sr4": decimal.Decimal("429228004229872.992467107"),
    "PTB-Yb1E3": decimal.Decimal("642121496772645.118522185"),
    "SYRTE-Hg": decimal.Decimal("1128575290808154.31910196"),
    # "SYRTE-Sr2": decimal.Decimal("429228004229872.992467107"),
    # "PTB-Yb1E3E2": decimal.Decimal("688358979309308.239120953"),
    "NPL-E3Yb+3": decimal.Decimal("642121496772645.118522185"),
    "NPL-Sr1": decimal.Decimal("429228004229872.992467107"),
}

print("TOCK March 2025")

# preprocessed link data to load
link_names = [
    # "INRIM_RioMod-SYRTE_CUS",
    # "NPL_T1-INRIM_RioMod",
    # "NPL_T1-SYRTE_CUS",
    # "PTB_NIRP-INRIM_RioMod",
    # "PTB_NIRP-NPL_T1",
    # "PTB_NIRP-SYRTE_CUS",
]

# other comparators to load
clock_names = [
    "INRIM_DoPTBSr4-INRIM_PTBSr4",
    "INRIM_DoPTBSr4-INRIM_LoYb",
    "INRIM_LoYb-INRIM_ITYb1",
]


start = 60740
stop = 60780


links = {}
for name in link_names:
    links[name] = rl.load_link_from_dir(
        os.path.join(outdir, name), start=ti.epoch_from_mjd(start), stop=ti.epoch_from_mjd(stop)
    )


for name in clock_names:
    # tstart, tstop = startstop[name]
    tstart, tstop = start, stop
    links[name] = rl.load_link_from_dir(
        os.path.join(dir, name), start=ti.epoch_from_mjd(tstart), stop=ti.epoch_from_mjd(tstop)
    )

# TODO
# I should make a list of the keys to fix the order!

chains = {
    "IT-Yb1/PTB-Sr4": ["INRIM_ITYb1", "INRIM_LoYb", "INRIM_DoPTBSr4", "INRIM_PTBSr4"],
}


# load all data
ratios = {}
for name, dos in chains.items():
    res = []
    for a, b in zip(dos[::-1][:-1], dos[::-1][1:]):
        if a + "-" + b in (link_names + clock_names):
            res += [-links[a + "-" + b]]
        elif b + "-" + a in (link_names + clock_names):
            res += [links[b + "-" + a]]
        else:
            print(a + "-" + b)

    ratios[name] = res


reslinks = {}

i = 0

for shortname, llinks in list(ratios.items()):
    print(shortname)
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


print(f"Saving clock data time = {time.time() - t0} s")
