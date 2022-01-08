#!/usr/bin/env python
# coding: utf-8

# In[95]:


import json
import logging

import holoviews as hv

from cdftpy import IonSolvation

from cdftpy.cdft1d.config import DATA_DIR
from cdftpy.cdft1d.rdf import analyze_rdf_peaks_sim
from cdftpy.cdft1d.solvent import Solvent

logging.getLogger('matplotlib').setLevel(logging.WARNING)

hv.extension('bokeh')
hv.opts.defaults(
    hv.opts.Layout( merge_tools=False, shared_axes=False, shared_datasource=False, framewise=True),
    hv.opts.Curve(width=500, height=500, show_grid=True, shared_axes=False, tools=['hover']),
    hv.opts.Scatter(width=500, height=500, show_grid=True, shared_axes=False, tools=['hover'])
)



with open(F"{DATA_DIR}/ion_data.json") as fp:
    ion_data = json.load(fp)
del ion_data["units"]
del ion_data["reference"]


with open(F"{DATA_DIR}/ion_md_data.json") as fp:
    ion_md_data = json.load(fp)


print(ion_md_data)


filename = F"{DATA_DIR}/s2.smdl"
solvent = Solvent.from_file(filename)

sim_array = []
names = []
charges = []
for name, pot in ion_data.items():
    charge = pot["charge"]
    if charge < 0:
        name = name + "-"
    elif charge >0:
        name = name + "+"
    charges.append(charge)
    names.append(name.capitalize())
    solute = {"name": name, **pot}
    _params = dict(ndiis=2, tol=1.0e-7, max_iter=500)
    sim = IonSolvation(**solute, solvent=solvent, method='rsdft',**_params)
    sim_array.append(sim)
    print(solute)
    fe = sim.cdft()
    print(fe)

    rdf_peaks = []
    fe_array = []
    for sim in sim_array:
        fe_array.append(sim.fe_tot)
        rdf_peaks.append(analyze_rdf_peaks_sim(sim))

o_peak_pos_md = []
o_peak_value_md = []
h_peak_pos_md = []
h_peak_value_md = []
o_peak_value = []
o_peak_pos = []
h_peak_value = []
h_peak_pos = []
for pk, name in zip(rdf_peaks,names):
    o_peak = pk["O"]['first_peak']
    o_peak_value.append(f"{round(o_peak[0],1)} ({ion_md_data[name]['O']['first_peak_value']})")
    o_peak_pos.append(f"{round(o_peak[1],1)} ({ion_md_data[name]['O']['first_peak_pos']})")
    h_peak = pk["H"]['first_peak']
    h_peak_value.append(f"{round(h_peak[0], 1)} ({ion_md_data[name]['H']['first_peak_value']})")
    h_peak_pos.append(f"{round(h_peak[1], 1)} ({ion_md_data[name]['H']['first_peak_pos']})")
    o_peak_value_md.append(ion_md_data[name]["O"]["first_peak_value"])
    o_peak_pos_md.append(ion_md_data[name]["O"]["first_peak_pos"])
    h_peak_value_md.append(ion_md_data[name]["H"]["first_peak_value"])
    h_peak_pos_md.append(ion_md_data[name]["H"]["first_peak_pos"])

table = hv.Table((names, fe , o_peak_value, o_peak_pos, h_peak_value, h_peak_pos),
                 ['Solute', 'Free Energy', 'Peak Height O ', 'Peak Pos O ','Peak Height H ', 'Peak Pos H ' ]
                ).opts(fontsize=14)
tbl = table.opts(height=340, width=500, title="Results of RSDFT simulations without the bridge")
hv.save(tbl, "table-clean.html")
