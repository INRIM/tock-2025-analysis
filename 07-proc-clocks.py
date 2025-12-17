import decimal
import os

import allantools as at
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.interpolate
import scipy.optimize
import tintervals as ti

# sys.path.append("..")
# import rocitlinks3 as rl
import tintervals.rocitlinks as rl
from matplotlib.backends.backend_pdf import PdfPages
from scipy.ndimage import minimum_filter1d

plt.close("all")
plt.ioff()

import time
from datetime import datetime

from scipy.ndimage import uniform_filter1d
from uncertainties import *

t0 = time.time()

dir = "./Data"
outdir = "./Analysis output"
figdir = os.path.join(outdir, "Figures")
stabdir = os.path.join(outdir, "Adev")
daydir = os.path.join(outdir, "Daily")

os.makedirs(outdir, exist_ok=True)
os.makedirs(figdir, exist_ok=True)
os.makedirs(stabdir, exist_ok=True)
os.makedirs(daydir, exist_ok=True)


def rlplot(link, *args, **kwargs):
    plt.plot(ti.mjd_from_epoch(link.t), link.delta, *args, **kwargs)


def rlplot_rolling(link, *args, **kwargs):
    plt.plot(ti.mjd_from_epoch(link.t), uniform_filter1d(link.delta, 1000), *args, **kwargs)


filo = open(os.path.join(outdir, "summary.dat"), "w")
filo.write("#Ratio		y(u)	y	uTot	uA*	uB1	uB2	uGRS	birge	ndof	days\n")

filoR = open(os.path.join(outdir, "ratios.dat"), "w")
filoR.write("# Ratio outputs\n")

pdfname = "INRIM_2025_analysis_" + datetime.today().strftime("%Y-%m-%d") + ".pdf"
pdf = PdfPages(pdfname)


print("TOCK March 2025")


names = np.loadtxt(os.path.join(outdir, "ratios.txt"), dtype=str)
ratio_names = {x[0]: x[1] for x in names}

# remote/local GRS uncertainty
clock_uGRS = {
    "IT-Yb1": (2.7, 1),
    "PTB-In1": (2.4, 0.33),
    "PTB-Al+": (2.4, 0.33),  # TBC
    "PTB-Yb1E3": (2.4, 0.31),
    "PTB-Sr3": (2.4, 0.39),
    "PTB-Sr4": (2.7, 1),
    "OBSPARIS-Sr2": (3, 1),
    "OBSPARIS-SrB": (3, 1),
    "PTB-Yb1E3E2": (2.4, 0.31),
    "NPL-E3Yb+3": (2.5, 0.38),
    "NPL-Sr1": (2.7, 1.14),
}


ratios = {}


start = 60842
stop = 60846

for shortname, longname in ratio_names.items():
    try:
        res = rl.load_link_from_dir(
            os.path.join(outdir, longname),
            discard_invalid=True,
            start=ti.epoch_from_mjd(start),
            stop=ti.epoch_from_mjd(stop),
        )
    except IOError:
        print(f"No data found for {shortname}")
        continue    

    if len(res.t) == 0:
        print(f"No data loaded for {shortname}")
        continue


    ratios[shortname] = res

extremes = np.array([[np.amin(x.t), np.amax(x.t)] for x in ratios.values()])
extremes = ti.mjd_from_epoch(extremes)

# link = ratios['IT-Yb1/PTB-Sr3']
# vals = [60768.567528, 60769.430905]
# mask = (link.t > ti.epoch_from_mjd(vals[0])) & (link.t < ti.epoch_from_mjd(vals[1]))
# link.flag[mask] = 0
# link.drop_invalid()

figs = []
i = 0

for shortname, reslink in list(ratios.items()):
    print(shortname)
    clock2, clock1 = shortname.split("/")
    print(reslink.name)
    print(reslink.r0)

    limit = (np.percentile(reslink.delta, 95) - np.percentile(reslink.delta, 5)) * 2
    print(f"{limit=:.2}")

    mask = np.abs(reslink.delta - np.mean(reslink.delta)) > limit
    reslink.flag[mask] = 0
    reslink.drop_invalid()
    outliers = np.sum(mask)

    if len(reslink.t) == 0:
        print("No common data!\n")
        continue

    if len(reslink.t) < 864:
        print("Less than 864 s of common data! Skipping!\n")
        continue

    print("Outliers = ", outliers)

    # time varying uB
    uuB1 = reslink.data[:, 3]
    uuB2 = reslink.data[:, 4]

    taua, ada, ade, adn = at.totdev(reslink.delta, rate=1, data_type="freq", taus="octave")
    tau, ad, ade, adn = at.oadev(reslink.delta, rate=1, data_type="freq", taus="octave")

    # calc uncertainty
    # https://tf.boulder.nist.gov/general/pdf/666.pdf
    # https://tf.nist.gov/general/pdf/666.pdf
    N = adn
    m = tau
    edf = (3 * (N - 1) / (2 * m) - 2 * (N - 2) / N) * 4 * m**2 / (4 * m**2 + 5)
    edf = np.clip(edf, 1, None)
    admax = ad * (edf / scipy.stats.chi2.ppf(0.1, edf)) ** 0.5
    admin = ad * (edf / scipy.stats.chi2.ppf(0.9, edf)) ** 0.5

    tottime = len(reslink.t) * reslink.step

    taufit = np.append(tau, tottime)

    # white calculated at 1000 s
    idx = (np.abs(tau - 1000)).argmin()
    white = ad[idx] * tau[idx] ** 0.5

    grsc = -reslink.oscA.grs_correction + reslink.oscB.grs_correction

    y = np.mean(reslink.delta) + grsc

    uA = white / tottime**0.5
    uuB = (uuB1**2 + uuB2**2) ** 0.5
    uB1 = np.mean(uuB1)
    uB2 = np.mean(uuB2)
    uB = (uB1**2 + uB2**2) ** 0.5

    fclock2, fclock1 = reslink.name.split("-")
    ist1 = fclock1.split("-")[0]
    ist2 = fclock2.split("-")[0]

    if ist1 == ist2:
        # Local GRS
        ug1, ug2 = clock_uGRS[clock1][1], clock_uGRS[clock1][1]
    else:
        ug1, ug2 = clock_uGRS[clock1][0], clock_uGRS[clock1][0]
    # ug1, ug2 = dict_uGRS[shortname]
    uGRS = (ug1**2 + ug2**2) ** 0.5 * 1e-18

    daily_vals = ti.array2intervals(reslink.t, tgap=3 * 3600)
    daily_vals = ti.semi_split(daily_vals, base=24 * 3600, offset=0 * 3600, mino=0.3)
    # dayly_vals = ti.array2intervals(reslink.t, tgap=3600)
    days, ddata, dcount = rl.average(reslink, daily_vals)
    mask = dcount > 864
    # mask = dcount > 86
    if sum(mask) > 0:
        days, ddata, dcount = days[mask], ddata[mask], dcount[mask]

    daily_ustat = np.sqrt(white**2 / (dcount))

    # resize fix the case that neither clock has variable uB given
    sys_vals, usys, sys_count = ti.maverage(np.resize(uuB, len(reslink.t)), reslink.t, days)

    daily_usys = (usys**2 + uGRS**2) ** 0.5

    ndof = len(daily_ustat) - 1
    if ndof > 0:
        chi2red = np.sum((ddata[:, 1] + grsc - y) ** 2 / daily_ustat**2) / ndof
        birge = np.sqrt(chi2red)
        uAstar = uA * max(birge, 1)
    else:
        birge = np.nan
        uAstar = uA

    final_u = np.sqrt(uAstar**2 + uB**2 + uGRS**2)
    final = ufloat(y, final_u)
    ratio = (decimal.Decimal(y) + decimal.Decimal(1)) * reslink.r0
    ratio_unc = float(ratio) * final_u

    # bug probne way to print the physical ratio
    def ratio_out(r, u):
        out = "{:0<18,.0f}({:.0f})".format(r * decimal.Decimal(1e18), u * 1e18)
        if r > 1:
            out = out[0] + "." + out[2:]
        else:
            out = "0." + out
        out = out.replace(",", " ")
        return out

    rout = ratio_out(ratio, ratio_unc)

    fig = plt.figure(figsize=(10, 6))
    gs = fig.add_gridspec(2, 2)

    ax0 = fig.add_subplot(gs[0, 0])
    fig.suptitle(
        f"{shortname}: Tot time: {tottime/3600:.2f} h, outliers = {outliers}\n y = {final:.2uS}  "
        + f"uA* = {uAstar:.1e}, uB1 = {uB1:.1e}, uB2 = {uB2:.1e},\n uGRS = {uGRS:.1e}, birge = {birge:.1f}  ndof = {ndof}\n R = "
        + rout
    )

    ax0.plot(ti.mjd_from_epoch(reslink.t) - 60000, reslink.delta + grsc, ".", label=reslink.name, rasterized=True)
    ax0.plot(
        ti.mjd_from_epoch(reslink.t) - 60000,
        uniform_filter1d(reslink.delta + grsc, 10000),
        ".",
        label="Rolling mean",
        rasterized=True,
    )
    ax0.legend(loc=0)
    ax0.set_ylabel("y")
    ax0.set_xlabel("MJD-60000")
    ax0.set_xlim(np.amin(extremes) - 60000, np.amax(extremes) - 60000)
    ax0.grid(True)

    # hours, hlink, hcount = rl.link_average(reslink, 360, timetags_as_start=False)
    # hlink.drop_invalid(0.5)

    # plt.figure()
    # # plt.errorbar(ti.mjd_from_epoch(hlink.t), hlink.delta, yerr=(white**2/(hlink.flag*hlink.step) + usys**2)**0.5, fmt='o', label=shortname)
    # # plt.errorbar(ti.mjd_from_epoch(hlink.t), hlink.delta, yerr=usys, fmt='o', label='Sys')
    # rlplot(hlink, '.', label=shortname)
    # plt.legend(loc=0)
    # plt.ylabel('y')
    # plt.xlabel('MJD')
    # plt.savefig(os.path.join(figdir, reslink.name + '.png'))

    # days, dlink, dcount = rl.link_average(reslink, 864, offset=0, timetags_as_start=False)
    # ustat = np.sqrt(white**2/(dlink.flag*dlink.step))

    # OPTION 1
    # dlink.drop_invalid(0.1)
    # timetags = ti.mjd_from_epoch(dlink.t)

    # OPTION 2
    timetags = ti.mjd_from_epoch(np.mean(days, axis=1))

    ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)
    ax1.axhspan(y - final_u, y + final_u, color="C2", alpha=0.5)
    # ax1.axhspan(y - (uA**2 + uB**2) ** 0.5, y + (uA**2 + uB**2) ** 0.5, color="C0", alpha=0.5)

    ax1.errorbar(
        timetags - 60000,
        ddata[:, 1] + grsc,
        yerr=(daily_ustat**2 + daily_usys**2) ** 0.5,
        fmt="o",
        label="Stat + Sys unc.",
    )
    ax1.errorbar(timetags - 60000, ddata[:, 1] + grsc, yerr=daily_ustat, fmt="o", label="Stat. unc.")
    # ax1.errorbar(timetags[outlier_mask], ddata[outlier_mask, 1] + grsc, fmt="o", label="Out")
    # ax1.plot(ti.mjd_from_epoch(dlink.t), dlink.delta  + grsc,  '.')
    ax1.set_ylabel("y (averaged)")
    ax1.set_xlabel("MJD-60000")
    ax1.set_xlim(np.amin(extremes) - 60000, np.amax(extremes) - 60000)
    ax1.legend(loc=0)
    ax1.grid(True)

    ax2 = fig.add_subplot(gs[:, 1])
    ax2.loglog(taufit, white / taufit**0.5, "-", label=f"White noise = {white:.2e}")
    ax2.errorbar(tau, ad, yerr=[ad - admin, admax - ad], fmt="o", label=reslink.name)
    ax2.loglog(taua, ada, "-", label="Tot dev")

    ax2.legend(loc=0)
    ax2.set_ylabel("ADEV")
    ax2.set_xlabel("Tau /s")
    ax2.grid(True, which="both")

    fig.tight_layout()

    # output
    print("Total time {:.2f} d".format(tottime / 86400))
    print("Mean {:.2uS}".format(final))

    print("Ratio " + rout)

    print("White noise at 1 s = {:.2}".format(white))
    print("Birge (white only) {:.2}".format(birge))
    print("utot= {:.2}".format(final_u))
    print("uA= {:.2}".format(uA))
    print("uA*= {:.2}".format(uAstar))
    print("uB1= {:.2}".format(uB2))  # note swap for consistency with ROCIT report
    print("uB2= {:.2}".format(uB1))
    print("uGRS= {:.2}".format(uGRS))
    print("\n")

    def round_sig(x, sig=2):
        if x != 0:
            return round(x, sig - int(np.floor(np.log10(abs(x)))) - 1)
        else:
            return x

    unit = 1
    filo.write(
        f"{shortname:16}	{final/unit:.2uS}	{y/unit:.3e}	{round_sig(final_u/unit):g}	{round_sig(uAstar/unit):g}	{round_sig(uB1/unit):g}	{round_sig(uB2/unit):g}	{round_sig(uGRS/unit):g}	{birge:.2f}	{ndof}	{tottime/86400:.1f}\n"
    )
    filoR.write(f"{shortname:16} = " + rout + "\n")

    # filoE.write(f'{shortname:16}\n')
    # filoE.write(ti.iso_from_epoch(reslink.t[0]) + '\t' + ti.iso_from_epoch(reslink.t[-1]) + '\n')
    # filoE.write(rout + f'\t{len(reslink.t)}\n')
    # filoE.write(f'{uA:.2e}\t{dark_hi:.2e}\n')
    # filoE.write(f'{(uB1**2 + uB2**2)**0.5:.2e}\t{dark_lo:.2e}\n\n')

    plt.savefig(os.path.join(figdir, reslink.name + "all.png"))
    pdf.savefig()
    i += 1

    figs += [fig]

    # save daily data
    out = np.column_stack(
        (ti.mjd_from_epoch(days), ti.mjd_from_epoch(ddata[:, 0]), ddata[:, 1] + grsc, daily_ustat, daily_usys)
    )
    np.savetxt(
        os.path.join(daydir, f"{i} " + reslink.name + ".dat"),
        out,
        fmt=["%.6f"] * 3 + ["%.6e"] + ["%.1e"] * 2,
        header="\t".join(["MJDstart", "MJDstop", "MJDbar", "y", "uStat", "uSys"]),
        delimiter="\t",
    )

    # save stability data for sharing
    stab_out = np.column_stack((tau, ad, admin, admax))
    np.savetxt(
        os.path.join(stabdir, f"{i} " + reslink.name + "_adev.dat"),
        stab_out,
        fmt=["%.1f", "%.3e", "%.3e", "%.3e"],
        delimiter="\t",
    )


# # calculate and save overlaps
# t_overlap = np.zeros((len(reslinks), len(reslinks)))

# # note that this preserves order, but it is not guranateed!
# ratio_list = list(reslinks.values())

# for i, i_ratio in enumerate(ratio_list):
#     for j, j_ratio in enumerate(ratio_list[: i + 1]):
#         common_t = np.intersect1d(i_ratio.t, j_ratio.t)
#         t_overlap[i, j] = len(common_t)
#         t_overlap[j, i] = len(common_t)

# np.savetxt(os.path.join(outdir, "overlap_matrix.dat"), t_overlap)


# plt.show(block=True)
filo.close()
filoR.close()
pdf.close()
print("Analysis time = ", time.time() - t0)
