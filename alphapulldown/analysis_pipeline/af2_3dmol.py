#
# Author: Grzegorz Chojnowski @ EMBL-Hamburg
#

import os, sys, re
import glob
import iotbx
import iotbx.pdb

import py3Dmol
from scitbx import matrix
from scitbx.math import superpose

def parse_pdbstring(pdb_string):

    # there may be issues with repeated BREAK lines, that we do not use here anyway
    pdb_string_lines = []
    for _line in pdb_string.splitlines():
        if re.match(r"^BREAK$", _line):
            continue
        pdb_string_lines.append(_line)

    # arghhhh, why do you guys keep changing the interface?
    inp = iotbx.pdb.input(source_info=None, lines=pdb_string_lines)
    try:
        return inp.construct_hierarchy(sort_atoms=False), inp.crystal_symmetry()
    except:
        return inp.construct_hierarchy(), inp.crystal_symmetry()


# -------------------------------------------------------------------------


def read_ph(ifname, selstr=None, verbose=True):

    if verbose:
        print(" ==> Parsing a PDB/mmCIF file: %s" % ifname)

    with open(ifname, "r") as ifile:
        ph, symm = parse_pdbstring(ifile.read())

    if selstr is None:
        return ph, symm

    sel_cache = ph.atom_selection_cache()
    isel = sel_cache.iselection

    phsel = ph.select(isel(selstr))

    return phsel, symm

def parse_results(output, color='lDDT', models=5,window_size=400,top=5):

    color_set=["lDDT", "rainbow", "chain"]
    if color not in color_set:
        print(f"Unknown color: {color}, will set to 'lDDT'.")
        color='lDDT'

    datadir = os.path.expanduser(output)

    pdb_fnames = sorted(glob.glob("%s/ranked*.pdb" % datadir))
    #print(pdb_fnames)
    # expected pdb filename: {parent}/<job>/ranked_{rank}.pdb

    assert top <= len(pdb_fnames)

    pdb_to_vis=[f'{output}/ranked_{rank}.pdb' for rank in range(top)]

    ph_array = []
    for idx, fn in enumerate(pdb_to_vis):
        _ph, _symm = read_ph(fn)
        if len(ph_array) > 0:
            _s = superpose.least_squares_fit(
                ph_array[0].atoms().extract_xyz(),
                _ph.atoms().extract_xyz(),
                method=["kearsley", "kabsch"][0],
            )
            rtmx = matrix.rt((_s.r, _s.t))
            _ph.atoms().set_xyz(new_xyz=rtmx * _ph.atoms().extract_xyz())

        ph_array.append(_ph)

    chain_ids = [_.id for _ in _ph.chains()]
    frames = min(models, len(ph_array))
    if len(chain_ids) > 1 :
        view = py3Dmol.view(
            js="https://3dmol.org/build/3Dmol.js",
            width=window_size,
            height=window_size * frames,
            viewergrid=(frames,1),
        )

        for idx, _ph in enumerate(ph_array):
            viewer = (idx,0)
            view.addModel(_ph.as_pdb_string(), "pdb", viewer=viewer)
            view.zoomTo(viewer=viewer)
            set_3dmol_styles(view, viewer, chain_ids=chain_ids, color=color)

        view.show()



def set_3dmol_styles(
    view,
    viewer,
    chain_ids=1,
    color=["lDDT", "rainbow", "chain"][0],
    show_sidechains=False,
    show_mainchains=False,
):

    """
    borrowed from colabfolds notebook at
    https://github.com/sokrypton/ColabFold/blob/main/colabfold/pdb.py
    """

    if color == "lDDT":
        view.setStyle(
            {
                "cartoon": {
                    "colorscheme": {
                        "prop": "b",
                        "gradient": "roygb",
                        "min": 0,
                        "max": 100,
                    }
                }
            },
            viewer=viewer,
        )
    elif color == "rainbow":
        view.setStyle({"cartoon": {"color": "spectrum"}}, viewer=viewer)

    elif color == "chain":
        for cid, color in zip(
            chain_ids,
            [
                "lime",
                "cyan",
                "magenta",
                "yellow",
                "salmon",
                "white",
                "blue",
                "orange",
                "black",
                "green",
                "gray",
            ]
            * 2,
        ):
            view.setStyle({"chain": cid}, {"cartoon": {"color": color}}, viewer=viewer)
    if show_sidechains:
        BB = ["C", "O", "N"]
        view.addStyle(
            {
                "and": [
                    {"resn": ["GLY", "PRO"], "invert": True},
                    {"atom": BB, "invert": True},
                ]
            },
            {"stick": {"colorscheme": f"WhiteCarbon", "radius": 0.3}},
        )
        view.addStyle(
            {"and": [{"resn": "GLY"}, {"atom": "CA"}]},
            {"sphere": {"colorscheme": f"WhiteCarbon", "radius": 0.3}},
        )
        view.addStyle(
            {"and": [{"resn": "PRO"}, {"atom": ["C", "O"], "invert": True}]},
            {"stick": {"colorscheme": f"WhiteCarbon", "radius": 0.3}},
        )
    if show_mainchains:
        BB = ["C", "O", "N", "CA"]
        view.addStyle(
            {"atom": BB},
            {"stick": {"colorscheme": f"WhiteCarbon", "radius": 0.3}},
            viewer=viewer,
        )
