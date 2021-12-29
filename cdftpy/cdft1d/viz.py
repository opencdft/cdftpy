#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Visualization Routines
"""

import io
import logging
from contextlib import redirect_stdout

import holoviews as hv
import numpy as np
import panel as pn
from prettytable import PrettyTable, PLAIN_COLUMNS

from cdftpy.cdft1d.rdf import analyze_rdf_peaks_sim
from cdftpy.cdft1d.rism import rism_1d
from cdftpy.cdft1d.rsdft import rsdft_1d
from cdftpy.cdft1d.solvent import solvent_model_locate, Solvent

logging.getLogger('matplotlib').setLevel(logging.WARNING)

hv.extension('bokeh')

css = '''
.bk.bk-tab {
  font-size: 120%;
}
'''

pn.extension(sizing_mode="stretch_both", raw_css=[css])

hv.opts.defaults(
    hv.opts.Layout(merge_tools=False, shared_axes=False, shared_datasource=False, framewise=True),
    hv.opts.Curve(show_grid=True, shared_axes=False, tools=['hover'], muted_alpha=0.,
                  fontsize={'title': 14, 'labels': 14, 'xticks': 10, 'yticks': 10}),
    hv.opts.Scatter(show_grid=True, muted_alpha=0., shared_axes=False, tools=['hover']))

water = ['red', 'blue', 'blue']


def epot_dashboard(sim):
    epot = np.sum(sim.epot_r, axis=0)
    qv = sim.solvent.charge
    rho_0 = sim.solvent.density
    rho = rho_0 * np.einsum("a,an->n", qv, sim.h_r + 1)

    r = sim.ifft.rgrid

    epot_plot = []
    xlim = r[-1] / 2
    epot_plot.append(hv.Curve((r, epot), 'r', 'Potential (au)',
                              group="Total Electric Potential").opts(xlim=(0, xlim), width=500))
    epot_plot.append(hv.Curve((r, rho), 'r', 'Density',
                              group="Total Charge Density").opts(xlim=(0, xlim), width=500))

    epot_widget = pn.Column(hv.Layout(epot_plot).cols(1), sizing_mode='stretch_both')
    return epot_widget


def multi_solute_rdf_dashboard(var, values, sim):
    aname_array = sim[0].solvent.aname

    rdf_widget = pn.Column()
    for i, aname in enumerate(aname_array):
        xlim = 0
        rdf_plot = []
        for j, s in enumerate(sim):
            pk = analyze_rdf_peaks_sim(s)
            xlim = int(max(xlim, pk[aname]["second_peak"][1] * 2))
            rdf = s.h_r[i, :] + 1
            r = s.ifft.rgrid
            curve = hv.Curve((r, rdf[:]), 'r', '\u03C1/\u03C1₀',
                             group=F"Solvent Density for {aname}", label=F"{var}={values[j]}")
            curve.opts(fontsize={'legend': 6})
            rdf_plot.append(curve)
        overlay = hv.Overlay(rdf_plot).opts(xlim=(0, xlim))
        overlay = overlay.options({'Curve': {'color': hv.Cycle('Category20b')}})
        rdf_widget.append(pn.panel(overlay))
    return rdf_widget


def multi_solute_energy_dashboard(var, values, sim):
    e_widget = pn.Column()

    e = []
    for s in sim:
        e.append(s.fe_tot)
        r = s.ifft.rgrid
    eref = e[0]
    ediff = [(x - eref) / 4.184 for x in e]
    e_plot = hv.Curve((values, ediff), F"{var}", 'Free Energy (kcal/mol)', group=F"Free Energy")
    e_plot = e_plot * hv.Scatter((values, ediff), F"{var}", 'Free Energy', group=F"Free Energy").opts(size=10,
                                                                                                      marker='circle')
    e_plot.opts(padding=0.1)
    e_widget.append(pn.panel(e_plot))

    tbl = PrettyTable()

    tbl.set_style(PLAIN_COLUMNS)

    tbl.add_column(var.capitalize(), values)
    tbl.add_column("Free Energy(kj/mol)", e)
    tbl.add_column("Free Energy Diff(kcal/mol)", ediff)

    tbl.align = "r"
    tbl.float_format = ".3"
    txt = tbl.get_string()
    e_table = pn.pane.HTML(F"<pre>{txt}</pre>")

    e_widget.append(pn.Column(e_table, scroll=True))
    return e_widget


def rdf_dashboard(sim):
    rdf = sim.h_r + 1
    r = sim.ifft.rgrid

    pk = analyze_rdf_peaks_sim(sim)

    rdf_plot = []
    xlim = 0
    for i, name in enumerate(sim.solvent.aname):
        rdf_plot.append(hv.Curve((r, rdf[i, :]), 'r', '\u03C1/\u03C1₀', group="Solvent Density", label=F" {name}"))
        xlim = int(max(xlim, pk[name]["second_peak"][1] * 4))

    overlay = hv.Overlay(rdf_plot)
    overlay = overlay.opts(xlim=(0, xlim))
    overlay = overlay.options({'Curve': {'color': hv.Cycle(water)}})
    rdf_panel = pn.panel(overlay)
    peaks = rdf_peaks_dashboard(sim)
    rdf_widget = pn.Column(rdf_panel, peaks, sizing_mode='stretch_both')

    return rdf_widget


def rdf_peaks_dashboard(sim):
    pk = analyze_rdf_peaks_sim(sim)

    tbl = PrettyTable()
    tbl.set_style(PLAIN_COLUMNS)

    tbl.field_names = ["Site", "1st peak pos/height",
                       "2nd peak pos/height", "1st min pos/height"]

    for i, name in enumerate(sim.solvent.aname):
        height1, pos1 = pk[name]["first_peak"]
        height2, pos2 = pk[name]["second_peak"]
        height3, pos3 = pk[name]["first_min"]
        tbl.add_row([name, F"{pos1.round(2)}/{height1.round(3)}",
                     F"{pos2.round(2)}/{height2.round(3)}",
                    F"{pos3.round(2)}/{height3.round(3)}"]
                    )
    tbl.align = "r"

    txt = tbl.get_string()
    return pn.pane.HTML(F"<pre>{txt}</pre>", height=100)


def xi_dashboard(sim):
    xi = sim.xi_r
    f = sim.f_r
    r = sim.ifft.rgrid

    aname = sim.solvent.aname
    pk = analyze_rdf_peaks_sim(sim)
    charge = sim.solvent.charge

    xi_plot = []
    f_plot = []
    xlim = 0
    for i, name in enumerate(sim.solvent.aname):
        ii = i - 1
        if charge[i] < 0:
            xlim = 6
        else:
            xlim = r[-1]
        xi_plot.append(
            hv.Curve((r, xi[i, :]), 'r', 'xi', group="Correlation Hole", label=F" {name}").opts(xlim=(0, xlim),
                                                                                                width=400))

    for i, name in enumerate(sim.solvent.aname):
        ii = i - 1
        if charge[ii] > 0:
            xlim = 6
        else:
            xlim = r[-1]
        f_plot.append(
            hv.Curve((r, f[ii, :]), 'r', 'f(r)', group="Mayer Function", label=F" {aname[ii]}").opts(xlim=(0, xlim),
                                                                                                     width=400))

    # hack to fix color assignment
    xi_plot = [x for _, x in sorted(zip(sim.solvent.charge, xi_plot), reverse=True)]
    f_plot = [x for _, x in sorted(zip(sim.solvent.charge, f_plot), reverse=True)]
    xi_widget = pn.Column(hv.Layout(xi_plot + f_plot).cols(2), sizing_mode='stretch_both')

    return xi_widget


def pmf_dashboard(sim):
    beta = sim.beta
    rdf = sim.h_r + 1

    r = sim.ifft.rgrid

    pk = analyze_rdf_peaks_sim(sim)

    pmf_plot = []
    xlim = 0

    for i, name in enumerate(sim.solvent.aname):
        ref_zero = np.min(rdf[i, :])
        rdf_adj = rdf[i, :] - ref_zero
        i0 = np.where(rdf_adj  > 1.e-4)[0][0]
        pmf = -np.log(rdf_adj[i0:] ) / beta
        pmf_plot.append(hv.Curve((r[i0:], pmf), 'r', 'PMF (kJ/mol)',
                                 group="Solvent-solute PMF",
                                 label=F" {name}"))
        xlim = int(max(xlim, pk[name]["second_peak"][1] * 2))

    overlay = hv.Overlay(pmf_plot)
    overlay = overlay.opts(xlim=(0, xlim))
    overlay = overlay.options({'Curve': {'color': hv.Cycle(water)}})
    pmf_widget = pn.Column(overlay, sizing_mode='stretch_both')

    return pmf_widget


def results_dashboard(sim):
    fe_tot = sim.fe_tot.round(3)
    solute = sim.solute
    solute_txt = F"{solute.name} charge={solute.charge} " \
                 F"sigma={solute.sigma} Å epsilon={solute.eps} kj/mol"
    solvent = sim.solvent
    f = io.StringIO()
    with redirect_stdout(f):
        solvent.report()
    s = f.getvalue()

    sim_params = F"Method: {sim.method}\n" \
                 F"Box size: {sim.rmax}\n" \
                 F"Temp: {solvent.temp}"
    txt = F"""
    <h2>Free energy of solvation:</h2> <pre>{fe_tot} kj/mol</pre>
    <hr>
    <h2>Simulation parameters:</h2><pre>{sim_params}</pre>
    <hr>
    <h2>Solute</h2><pre>{solute_txt} </pre>
    <hr>
    <h2>Solvent</h2><pre>{s} </pre>
    """

    html_pane = pn.pane.HTML(F"{txt}")

    return pn.Column(html_pane, scroll=True, width=500)


def multi_solute_results_dashboard(sim):
    s = sim[0]
    solute = s.solute
    solute_txt = F"{solute.name} charge={solute.charge} " \
                 F"sigma={solute.sigma} Å epsilon={solute.eps} kj/mol"
    solvent = s.solvent
    f = io.StringIO()
    with redirect_stdout(f):
        solvent.report()
    solv_txt = f.getvalue()

    sim_params = F"Method: {s.method}\n" \
                 F"Box size: {s.rmax}\n" \
                 F"Temp: {solvent.temp}"
    txt = F"""
    <h2>Simulation parameters:</h2><pre>{sim_params}</pre>
    <hr>
    <h2>Solute</h2><pre>{solute_txt} </pre>
    <hr>
    <h2>Solvent</h2><pre>{solv_txt} </pre>
    """
    html_pane = pn.pane.HTML(F"{txt}")
    return pn.Column(html_pane, scroll=True, width=500)


def single_point_viz(sim, dashboard_dest="browser"):
    html_pane = results_dashboard(sim)
    rdf_widget = rdf_dashboard(sim)
    pmf_widget = pmf_dashboard(sim)
    epot_widget = epot_dashboard(sim)
    charts = pn.Tabs(("Density", rdf_widget), ("PMF", pmf_widget), ("Electric Potential", epot_widget))
    viz = pn.Row(html_pane, charts)
    if dashboard_dest == "browser":
        viz.show()
    else:
        viz.save(dashboard_dest)


def multi_solute_viz(var, values, sim, dashboard_dest="browser"):
    html_pane = multi_solute_results_dashboard(sim)
    rdf_widget = multi_solute_rdf_dashboard(var, values, sim)
    e_widget = multi_solute_energy_dashboard(var, values, sim)
    charts = pn.Tabs(("Solvation Free Energy", e_widget), ("Density", rdf_widget))

    viz = pn.Row(html_pane, charts)
    if dashboard_dest == "browser":
        viz.show()
    else:
        viz.save(dashboard_dest)


def single_solute_test():

    params = dict(diis_iterations=2, tol=1.0e-7, max_iter=200)
    solvent_name = "spce"
    filename = solvent_model_locate(solvent_name)
    solv0 = Solvent.from_file(filename, rism_patch=False)
    solute = dict(name="Cl", charge=-1.0, sigma=4.83, eps=0.05349244)
    solute = {
        "name": "F",
        "charge": -1.0,
        "sigma": 4.02,
        "eps": 0.030963692
    }
    # sim = rsdft_1d(solute, solv0, params=params, quiet=True)
    sim = rism_1d(solute, solv0, params=params, quiet=True)

    rdf_peaks_dashboard(sim)
    html_pane = results_dashboard(sim)
    rdf_widget = rdf_dashboard(sim)
    pmf_widget = pmf_dashboard(sim)
    epot_widget = epot_dashboard(sim)
    # xi_widget = xi_dashboard(sim)
    # charts = pn.Tabs(("Density",rdf_widget),("PMF",pmf_widget),("Correlation Hole",xi_widget))
    charts = pn.Tabs(("Density", rdf_widget), ("PMF", pmf_widget), ("Electric Potential", epot_widget))
    # pn.Row(html_pane,charts).save("analysis.html")
    pn.Row(html_pane, charts).show()
    #
    # rdf_widget.show()


def multi_solute_test():

    parameters = dict(diis_iterations=2, tol=1.0e-7, max_iter=200)
    solvent_name = "s2"
    filename = solvent_model_locate(solvent_name)
    solvent = Solvent.from_file(filename, rism_patch=False)
    solute = dict(name="Cl", charge=-1.0, sigma=4.83, eps=0.05349244)

    sim = []
    values =  [-1, -0.5, 0]
    for v in values:
        solute["charge"] = v
        s = rsdft_1d(solute, solvent, params=parameters)
        sim.append(s)

    multi_solute_viz("charge", values, sim, dashboard_dest="browser")

if __name__ == '__main__':
    single_solute_test()