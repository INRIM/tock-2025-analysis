# 05-calc-link
# Calc link data with rocitlinks
# Copyright (C) 2022  Marco Pizzocaro - Istituto Nazionale di Ricerca Metrologica (INRIM)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# SPDX-License-Identifier: MIT


import decimal
import os
import time

import numpy as np
import tintervals as ti
import tintervals.rocitlinks as rl

t0 = time.time()


syrte_dos = ["OBSPARIS_CUS_2", "OBSPARIS_MLS"]
join_dos = ["OBSPARIS_MLS", "TH2_MLS"]
inrim_dos = ["TH2_MLS", "MODANE_RLS", "INRIM_RioMod"]
ptb_dos = ["TH2_MLS", "UNISTRA_RLS2", "UNISTRA_RLS5", "PTB_NIRP"]
npl_dos = ["OBSPARIS_MLS", "LPL_MLS", "LPL_RLS", "NPL_MLS", "NPL_T1"]


dir = "./Data"

start = 60748.0
stop = 60780.0

syrte = rl.load_links_from_osc_names(syrte_dos, dir, start=ti.epoch_from_mjd(start), stop=ti.epoch_from_mjd(stop))
join = rl.load_links_from_osc_names(join_dos, dir, start=ti.epoch_from_mjd(start), stop=ti.epoch_from_mjd(stop))
inrim = rl.load_links_from_osc_names(inrim_dos, dir, start=ti.epoch_from_mjd(start), stop=ti.epoch_from_mjd(stop))
# npl = rl.load_links_from_osc_names(npl_dos, dir, start=ti.epoch_from_mjd(start), stop=ti.epoch_from_mjd(stop))
ptb = rl.load_links_from_osc_names(ptb_dos, dir, start=ti.epoch_from_mjd(start), stop=ti.epoch_from_mjd(stop))

# this is not exposed in the yaml files
# but syrte cus ans INRIM riomod are decent reference oscillator
syrte[0].oscA.v0 = decimal.Decimal(194.4e12)
inrim[-1].oscB.v0 = decimal.Decimal(194.4e12)

# let's hope it is true also for NPL
# npl[-1].oscB.v0 = decimal.Decimal(194.4e12)


syrte_inrim = syrte + join + inrim
syrte_ptb = syrte + join + ptb
# syrte_npl = syrte + npl
inrim_ptb = [-x for x in inrim[::-1]] + ptb
# inrim_npl = [-x for x in inrim[::-1]] + [-join[0]] + npl
# npl_ptb = [-x for x in npl[::-1]] + join + ptb


outdir = "./Analysis output"
os.makedirs(outdir, exist_ok=True)


list_links = [
    syrte_inrim,
    syrte_ptb,
    #    syrte_npl,
    inrim_ptb,
    #    inrim_npl,
    #    npl_ptb,
]

for i, links in enumerate(list_links):
    reslink, masks = rl.chain(*links)

    rl.save_link_to_dir(outdir, reslink)


print("Link calc time = ", time.time() - t0)
