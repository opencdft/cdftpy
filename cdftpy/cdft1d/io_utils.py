#!/usr/bin/python
# -*- coding: utf-8 -*-
import ast
import sys
import numpy as np

def print_banner(title):
    dash = "-" * len(title)
    print(dash)
    print(title)
    print(dash)

def get_banner(title):
    dash = "-" * len(title)
    return F"{dash}\n{title}\n{dash}"


def iter_lines(fp):
    for line in fp:
        record = line.rsplit("#")[0].strip()
        if record == "":
            continue
        else:
            yield record


def read_key_value(filename, section):
    tag = f"<{section}>"

    params = {}

    with open(filename, "r") as fp:
        if find_section(tag, fp) is None:
            return None
        for line in iter_lines(fp):
            if line.startswith("<"):
                break
            tokens = line.split()
            if len(tokens) == 2:
                key, value = tokens
                try:
                    value = ast.literal_eval(value)
                except ValueError:
                    pass
            elif len(tokens) == 1:
                key, value = tokens[0], True
            params[key] = value
        return params


def read_solute(filename, section="solute"):
    name_list = []
    sigma_list = []
    eps_list = []
    charge_list = []


    tag = f"<{section}>"
    with open(filename, "r") as fp:
        if find_section(tag, fp) is None:
            print(f"cannot find {tag} section")
            sys.exit(1)
        for line in iter_lines(fp):
            if line.startswith("<"):
                break
            name, sigma, eps, charge, *xyz = line.split()
            name_list.append(name)
            sigma_list.append(float(sigma))
            eps_list.append(float(eps))
            charge_list.append(float(charge))


    return dict(
        name=np.array(name_list),
        sigma=np.array(sigma_list),
        eps=np.array(eps_list),
        charge=np.array(charge_list)
    )


def read_geom(filename, section="geometry"):
    name_list = []
    atype_list = []
    x_list = []
    y_list = []
    z_list = []

    tag = f"<{section}>"
    with open(filename, "r") as fp:
        if find_section(tag, fp) is None:
            print(f"cannot find {tag} section")
            sys.exit(1)
        for line in iter_lines(fp):
            if line.startswith("<"):
                break
            name, atype, x, y, z = line.split()
            name_list.append(name)
            atype_list.append(atype)
            x_list.append(float(x))
            y_list.append(float(y))
            z_list.append(float(z))

    return dict(
        aname=np.array(name_list),
        atype=np.array(atype_list),
        x=np.array(x_list),
        y=np.array(y_list),
        z=np.array(z_list),
    )


def read_interaction(filename, section="interactions"):
    ff = dict()
    tag = f"<{section}>"
    with open(filename, "r") as fp:
        if find_section(tag, fp) is None:
            print(f"cannot find {tag} section")
            sys.exit(1)
        for line in iter_lines(fp):
            if line.startswith("<"):
                break
            name, sigma, eps, charge = line.split()
            ff[name] = dict(sigma=float(sigma), eps=float(eps), charge=float(charge))

    return ff


def read_rism_patch(filename, section="rism_patch"):
    ff = dict()
    tag = f"<{section}>"
    with open(filename, "r") as fp:
        if find_section(tag, fp) is None:
            # print(f"cannot find {tag} section")
            return ff
        for line in iter_lines(fp):
            if line.startswith("<"):
                break
            name, sigma, eps = line.split()
            ff[name] = dict(sigma=float(sigma), eps=float(eps))

    return ff


def read_molecule(
        filename,
        geometry_section="geometry",
        interaction_section="interactions",
        rism_patch=False,
):
    sigma_list = []
    eps_list = []
    charge_list = []

    system = read_geom(filename, section=geometry_section)
    ff = read_interaction(filename, interaction_section)

    for i, name in enumerate(system["atype"]):
        sigma_list.append(ff[name]["sigma"])
        eps_list.append(ff[name]["eps"])
        charge_list.append(ff[name]["charge"])

    if rism_patch:
        ff_rism = read_rism_patch(filename)
        for i, name in enumerate(system["atype"]):
            if name in ff_rism:
                sigma_list[i] = ff_rism[name]["sigma"]
                eps_list[i] = ff_rism[name]["eps"]

    system["sigma"] = np.array(sigma_list)
    system["eps"] = np.array(eps_list)
    system["charge"] = np.array(charge_list)

    return system


def read_array(filename, section):
    tag = f"<{section}>"

    with open(filename, "r") as fp:
        if find_section(tag, fp) is None:
            # print(f"cannot find {tag} section in {os.path.abspath(filename)}")
            return None, None
        line = fp.readline()
        nv, ngrid = map(int, line.split())
        kgrid = np.zeros(shape=ngrid)
        a_k = np.zeros(shape=(nv, nv, ngrid), dtype=np.double)

        for ig, line in enumerate(iter_lines(fp)):
            if line.startswith("<"):
                break
            row = list(map(float, line.split()))
            kgrid[ig] = row[0]
            n = 1
            for i in range(nv):
                for j in range(i, nv):
                    a_k[i, j, ig] = row[n]
                    a_k[j, i, ig] = row[n]
                    n = n + 1

        return a_k, kgrid


def find_section(section, fp):
    iline = None
    for i, line in enumerate(fp, 1):
        if line.strip() == section:
            iline = i
            break

    return iline


# def sim_to_netcdf(sim):
#
#     rdims = ["site", "rgrid"]
#     kdims = ["site", "kgrid"]
#     data_r = dict(h_r=(rdims, sim.h_r),
#                   g_r=(rdims, sim.g_r),
#                   gl_r=(rdims, sim.gl_r),
#                   vs_r=(rdims, sim.vs_r),
#                   cs_r=(rdims, sim.cs_r),
#                   xi_r=(rdims, sim.xi_r),
#                   epot_r=(rdims, sim.epot_r),
#                   vl_r=(rdims, sim.vl_r)
#                   )
#
#     data_k = dict(h_k=(kdims, sim.h_k),
#                   gl_k=(kdims, sim.gl_k),
#                   vl_k=(kdims, sim.vl_k),
#                   cs_k=(kdims, sim.cs_k),
#                   epot_k=(kdims, sim.epot_k)
#                   )
#
#     coords = dict(site=sim.solvent.aname, rgrid=sim.rgrid, kgrid=sim.kgrid)
#     ds = xr.Dataset(data_vars=data_r | data_k, coords=coords, attrs=dict(solute=solute))
#     print(ds)
#     ds.to_netcdf("test.nc", engine='h5netcdf')

