"""
netParams.py 

High-level specifications for M1 network model using NetPyNE

Contributors: salvadordura@gmail.com
"""

from netpyne import specs
import pickle, json
import glob
import copy
import numpy as np
from netpyne.specs.dicts import Dict

netParams = (
    specs.NetParams()
)  # object of class NetParams to store the network parameters

netParams.version = 56

try:
    from __main__ import cfg  # import SimConfig object with params from parent module
except:
    from cfg import cfg


# Take from BluePyOpt
def replace_axon(cell, axon_stub_length=60, axon_nseg_frequency=40):
    """Replace axon"""

    axon = [cell.secs[sec] for sec in cell.secs if "axon" in sec]
    nsec = len(axon)
    if "soma" in cell.secs:
        soma = cell.secs["soma"]
    elif "soma_0":
        soma = cell.secs["soma_0"]
    else:
        raise Exception("cell missing soma")
    if nsec == 0:
        ais_diams = [1, 1]
    elif nsec == 1:
        ais_diams = [axon[0]["geom"]["diam"], axon[0]["geom"]["diam"]]
    else:
        ais_diams = [axon[0]["geom"]["diam"], axon[0]["geom"]["diam"]]

        dist = cell["soma"]["geom"]["L"] / 2
        for sec in axon:
            dist += sec["geom"]["L"]
            if axon_stub_length:
                ais_diams[1] = sec["geom"]["diam"]
                break

    # Only keep 1 section
    remove = False
    for sec in cell.secs.keys():
        if "axon" in sec:
            if remove:
                del cell.secs[sec]
            else:
                remove = True

    # new geometry
    L = axon_stub_length
    nseg = 1 + 2 * int(L / axon_nseg_frequency)
    if nseg == 1:
        diams = ais_diams
    else:
        x = np.linspace(0, L, nseg)
        diams = np.interp(x, [0, L], ais_diams)

    # calculate new pts3d
    pts = np.array(axon[0]["geom"]["pt3d"])
    x, y, z = pts[0, :3] - pts[-1, :3]
    Lorig = sum((pts[0, :3] - pts[-1, :3]) ** 2) ** 0.5
    x0, y0, z0, _ = soma["geom"]["pt3d"][-1]
    newpts = [(x0, y0, z0, ais_diams[0])]
    dx, dy, dz = (
        L * x / Lorig / nseg,
        L * y / Lorig / nseg,
        L * z / Lorig / nseg,
    )
    pt = newpts[0]
    for i in range(nseg):
        a, b, c, _ = pt
        pt = (a + dx, b + dy, c + dz, diams[i])
        newpts.append(pt)
    axon[0]["geom"]["L"] = L
    axon[0]["geom"]["nseg"] = nseg
    axon[0]["geom"]["diam"] = ais_diams[0]  # if nseg == 1 else diams
    axon[0]["geom"]["pt3d"] = newpts
    axon[0]["topol"]["parentX"] = 1.0


# ------------------------------------------------------------------------------
#
# NETWORK PARAMETERS
#
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# General network parameters
# ------------------------------------------------------------------------------
netParams.scale = cfg.scale  # Scale factor for number of cells
netParams.sizeX = cfg.sizeX  # x-dimension (horizontal length) size in um
netParams.sizeY = (
    cfg.sizeY
)  # y-dimension (vertical height or cortical depth) size in um
netParams.sizeZ = cfg.sizeZ  # z-dimension (horizontal depth) size in um
netParams.shape = "cylinder"  # cylindrical (column-like) volume

# ------------------------------------------------------------------------------
# General connectivity parameters
# ------------------------------------------------------------------------------
netParams.scaleConnWeight = (
    1.0  # Connection weight scale factor (default if no model specified)
)
netParams.scaleConnWeightModels = {
    "HH_simple": 1.0,
    "HH_reduced": 1.0,
    "HH_full": 1.0,
}  # scale conn weight factor for each cell model
netParams.scaleConnWeightNetStims = 1.0  # 0.5  # scale conn weight factor for NetStims
netParams.defaultThreshold = (
    0.0  # spike threshold, 10 mV is NetCon default, lower it for all cells
)
netParams.defaultDelay = 2.0  # default conn delay (ms)
netParams.propVelocity = 500.0  # propagation velocity (um/ms)
netParams.probLambda = (
    100.0  # length constant (lambda) for connection probability decay (um)
)
netParams.defineCellShapes = True  # convert stylized geoms to 3d points


# special condition to change Kgbar together with ih when running batch
# note min Kgbar is assumed to be 0.5, so this is set here as an offset
if cfg.makeKgbarFactorEqualToNewFactor:
    cfg.KgbarFactor = 0.5 + cfg.modifyMechs["newFactor"]

# ------------------------------------------------------------------------------
# Cell parameters
# ------------------------------------------------------------------------------
cellModels = ["HH_simple", "HH_reduced", "HH_full"]
layer = {
    "1": [0.0, 0.1],
    "2": [0.1, 0.29],
    "4": [0.29, 0.37],
    "5A": [0.37, 0.47],
    "24": [0.1, 0.37],
    "5B": [0.47, 0.8],
    "6": [0.8, 1.0],
    "longTPO": [2.0, 2.1],
    "longTVL": [2.1, 2.2],
    "longS1": [2.2, 2.3],
    "longS2": [2.3, 2.4],
    "longcM1": [2.4, 2.5],
    "longM2": [2.5, 2.6],
    "longOC": [2.6, 2.7],
}  # normalized layer boundaries

netParams.correctBorder = {
    "threshold": [
        cfg.correctBorderThreshold,
        cfg.correctBorderThreshold,
        cfg.correctBorderThreshold,
    ],
    "yborders": [layer["2"][0], layer["5A"][0], layer["6"][0], layer["6"][1]],
}  # correct conn border effect

# ------------------------------------------------------------------------------
## Load cell rules previously saved using netpyne format
cellParamLabels = [
    "IT2_reduced",
    "IT4_reduced",
    "IT5A_reduced",
    "IT5B_reduced",
    "PT5B_reduced",
    "IT6_reduced",
    "CT6_reduced",
    "PV_simple",
    "SOM_simple",
    "IT5A_full",
]  # 'VIP_reduced', 'NGF_simple','PT5B_full'] #  # list of cell rules to load from file
loadCellParams = cellParamLabels
saveCellParams = False  # True

for ruleLabel in loadCellParams:
    netParams.loadCellParamsRule(
        label=ruleLabel, fileName="cells/" + ruleLabel + "_cellParams.pkl"
    )

    # Adapt K gbar
    if ruleLabel in [
        "IT2_reduced",
        "IT4_reduced",
        "IT5A_reduced",
        "IT5B_reduced",
        "IT6_reduced",
        "CT6_reduced",
        "IT5A_full",
    ]:
        cellRule = netParams.cellParams[ruleLabel]
        for secName in cellRule["secs"]:
            for kmech in [
                k
                for k in cellRule["secs"][secName]["mechs"].keys()
                if k.startswith("k") and k != "kBK"
            ]:
                cellRule["secs"][secName]["mechs"][kmech]["gbar"] *= cfg.KgbarFactor

# ------------------------------------------------------------------------------
# Specification of cell rules not previously loaded
# Includes importing from hoc template or python class, and setting additional params

# ------------------------------------------------------------------------------
# Reduced cell model params (6-comp)
reducedCells = {  # layer and cell type for reduced cell models
    "IT2_reduced": {"layer": "2", "cname": "CSTR6", "carg": "BS1578"},
    "IT4_reduced": {"layer": "4", "cname": "CSTR6", "carg": "BS1578"},
    "IT5A_reduced": {"layer": "5A", "cname": "CSTR6", "carg": "BS1579"},
    "IT5B_reduced": {"layer": "5B", "cname": "CSTR6", "carg": "BS1579"},
    "PT5B_reduced": {"layer": "5B", "cname": "SPI6", "carg": None},
    "IT6_reduced": {"layer": "6", "cname": "CSTR6", "carg": "BS1579"},
    "CT6_reduced": {"layer": "6", "cname": "CSTR6", "carg": "BS1578"},
}

reducedSecList = {  # section Lists for reduced cell model
    "alldend": ["Adend1", "Adend2", "Adend3", "Bdend"],
    "spiny": ["Adend1", "Adend2", "Adend3", "Bdend"],
    "apicdend": ["Adend1", "Adend2", "Adend3"],
    "perisom": ["soma"],
}

for label, p in reducedCells.items():  # create cell rules that were not loaded
    if label not in loadCellParams:
        cellRule = netParams.importCellParams(
            label=label,
            conds={
                "cellType": label[0:2],
                "cellModel": "HH_reduced",
                "ynorm": layer[p["layer"]],
            },
            fileName="cells/" + p["cname"] + ".py",
            cellName=p["cname"],
            cellArgs={"params": p["carg"]} if p["carg"] else None,
        )
        dendL = (
            layer[p["layer"]][0] + (layer[p["layer"]][1] - layer[p["layer"]][0]) / 2.0
        ) * cfg.sizeY  # adapt dend L based on layer
        for secName in ["Adend1", "Adend2", "Adend3", "Bdend"]:
            cellRule["secs"][secName]["geom"]["L"] = dendL / 3.0  # update dend L
        for k, v in reducedSecList.items():
            cellRule["secLists"][k] = v  # add secLists
        netParams.addCellParamsWeightNorm(
            label,
            "conn/" + label + "_weightNorm.pkl",
            threshold=cfg.weightNormThreshold,
        )  # add weightNorm

        # set 3d points
        offset, prevL = 0, 0
        somaL = netParams.cellParams[label]["secs"]["soma"]["geom"]["L"]
        for secName, sec in netParams.cellParams[label]["secs"].items():
            sec["geom"]["pt3d"] = []
            if secName in [
                "soma",
                "Adend1",
                "Adend2",
                "Adend3",
            ]:  # set 3d geom of soma and Adends
                sec["geom"]["pt3d"].append([offset + 0, prevL, 0, sec["geom"]["diam"]])
                prevL = float(prevL + sec["geom"]["L"])
                sec["geom"]["pt3d"].append([offset + 0, prevL, 0, sec["geom"]["diam"]])
            if secName in ["Bdend"]:  # set 3d geom of Bdend
                sec["geom"]["pt3d"].append([offset + 0, somaL, 0, sec["geom"]["diam"]])
                sec["geom"]["pt3d"].append(
                    [offset + sec["geom"]["L"], somaL, 0, sec["geom"]["diam"]]
                )
            if secName in ["axon"]:  # set 3d geom of axon
                sec["geom"]["pt3d"].append([offset + 0, 0, 0, sec["geom"]["diam"]])
                sec["geom"]["pt3d"].append(
                    [offset + 0, -sec["geom"]["L"], 0, sec["geom"]["diam"]]
                )

        if saveCellParams:
            netParams.saveCellParamsRule(
                label=label, fileName="cells/" + label + "_cellParams.pkl"
            )

# ------------------------------------------------------------------------------
## PT5B model
if cfg.control is None:  # load the original model
    ## PT5B full cell model params (700+ comps)
    if "PT5B_full" not in loadCellParams:
        ihMod2str = {"harnett": 1, "kole": 2, "migliore": 3}
        cellRule = netParams.importCellParams(
            label="PT5B_full",
            conds={"cellType": "PT", "cellModel": "HH_full"},
            fileName="cells/PTcell.hoc",
            cellName="PTcell",
            cellArgs=[ihMod2str[cfg.ihModel], cfg.ihSlope],
            somaAtOrigin=True,
        )
        nonSpiny = ["apic_0", "apic_1"]
        netParams.addCellParamsSecList(
            label="PT5B_full", secListName="perisom", somaDist=[0, 50]
        )  # sections within 50 um of soma
        netParams.addCellParamsSecList(
            label="PT5B_full", secListName="below_soma", somaDistY=[-600, 0]
        )  # sections within 0-300 um of soma
        for sec in nonSpiny:
            cellRule["secLists"]["perisom"].remove(sec)
        cellRule["secLists"]["alldend"] = [
            sec for sec in cellRule.secs if ("dend" in sec or "apic" in sec)
        ]  # basal+apical
        cellRule["secLists"]["apicdend"] = [
            sec for sec in cellRule.secs if ("apic" in sec)
        ]  # apical
        cellRule["secLists"]["spiny"] = [
            sec for sec in cellRule["secLists"]["alldend"] if sec not in nonSpiny
        ]
        # Adapt ih params based on cfg param
        for secName in cellRule["secs"]:
            for mechName, mech in cellRule["secs"][secName]["mechs"].items():
                if mechName in ["ih", "h", "h15", "hd"]:
                    mech["gbar"] = (
                        [g * cfg.ihGbar for g in mech["gbar"]]
                        if isinstance(mech["gbar"], list)
                        else mech["gbar"] * cfg.ihGbar
                    )
                    if cfg.ihModel == "migliore":
                        mech["clk"] = cfg.ihlkc  # migliore's shunt current factor
                        mech["elk"] = (
                            cfg.ihlke
                        )  # migliore's shunt current reversal potential
                    if secName.startswith("dend"):
                        mech[
                            "gbar"
                        ] *= (
                            cfg.ihGbarBasal
                        )  # modify ih conductance in soma+basal dendrites
                        mech[
                            "clk"
                        ] *= (
                            cfg.ihlkcBasal
                        )  # modify ih conductance in soma+basal dendrites
                    if (
                        secName in cellRule["secLists"]["below_soma"]
                    ):  # secName.startswith('dend'):
                        mech[
                            "clk"
                        ] *= (
                            cfg.ihlkcBelowSoma
                        )  # modify ih conductance in soma+basal dendrites

            # Adapt K gbar
            for kmech in [
                k
                for k in cellRule["secs"][secName]["mechs"].keys()
                if k.startswith("k") and k != "kBK"
            ]:
                cellRule["secs"][secName]["mechs"][kmech]["gbar"] *= cfg.KgbarFactor

        ##### Adding params to change in batch ###################################
        #    cellRule['secs']['soma']['mechs']['nap']['gbar'] = cfg.gbar        # change gbar in nap channel (soma only)
        cellRule["secs"]["soma"]["mechs"]["kBK"][
            "gpeak"
        ] = cfg.gpeak  # change gpeak in kBK channel (soma only)
        #    cellRule['secs']['soma']['mechs']['Nafx']['gnafbar'] = cfg.gnafbar # change gnafbar in NaT channel (soma only)
        for secName in cellRule["secs"]:
            cellRule["secs"][secName]["mechs"]["nap"][
                "gbar"
            ] = cfg.gbar  # change gbar in nap channel (all sections)
        for secName in cellRule["secs"]:
            cellRule["secs"][secName]["mechs"]["Nafx"][
                "gnafbar"
            ] = cfg.gnafbar  # change gnafbar in NaT channel (all sections)

        # Reduce dend Na to avoid dend spikes (compensate properties by modifying axon params)
        for secName in cellRule["secLists"]["alldend"]:
            cellRule["secs"][secName]["mechs"]["nax"]["gbar"] = (
                0.0153130368342 * cfg.dendNa
            )  # 0.25

        cellRule["secs"]["soma"]["mechs"]["nax"]["gbar"] = 0.0153130368342 * cfg.somaNa
        cellRule["secs"]["axon"]["mechs"]["nax"]["gbar"] = (
            0.0153130368342 * cfg.axonNa
        )  # 11
        cellRule["secs"]["axon"]["geom"]["Ra"] = 137.494564931 * cfg.axonRa  # 0.005
        # Remove Na (TTX)
        if cfg.removeNa:
            for secName in cellRule["secs"]:
                cellRule["secs"][secName]["mechs"]["nax"]["gbar"] = 0.0
        netParams.addCellParamsWeightNorm(
            "PT5B_full",
            "conn/PT5B_full_weightNorm.pkl",
            threshold=cfg.weightNormThreshold,
        )  # load weight norm
        if saveCellParams:
            netParams.saveCellParamsRule(
                label="PT5B_full", fileName="cells/PT5B_full_cellParams.pkl"
            )
else:
    ## PT5B reduced model fit bit BluePyOpt
    PT5B = netParams.importCellParams(
        label="PT5B",
        conds={"cellType": "PT", "cellModel": "HH_reduced"},
        fileName="./cells/PT5B/HHLong.swc",
        cellName="PT5B_BPO",
    )
    replace_axon(PT5B)
    # add all mechanims
    mechs = json.load(open("./cells/PT5B/mechanisms_alt.json", "r"))
    for sec in PT5B["secs"]:
        sl = sec.split("_")[0]
        PT5B["secs"][sec]["mechs"] = {mname: {} for mname in mechs["all"]}
        for mname in mechs[sl]:
            PT5B["secs"][sec]["mechs"][mname] = {}
    netParams.renameCellParamsSec(
        "PT5B", "soma_0", "soma"
    )  # rename imported section 'soma_0' to 'soma'

    vinit = json.load(open("./cells/PT5B/vinit.json", "r"))

    # set the Nernst potentials
    for sec in PT5B["secs"]:
        PT5B["secs"][sec]["ions"]["na"]["e"] = cfg.pt5bEna
        PT5B["secs"][sec]["ions"]["k"]["e"] = cfg.pt5bEk

    # Add sectlists
    # spiny is used for synaptic connections
    reducedSecList = {
        "spiny": ["apic_0", "dend_0"],
        "perisom": ["soma"],
        "apicdend": ["apic_0"],
        "alldend": ["dend_0", "apic_0"],
    }
    for k, v in reducedSecList.items():
        PT5B["secLists"][k] = v

    def paramsBySec(params):
        pbysec = {}
        geom = {}
        for p, v in params.items():
            param, sec = p.split(".")
            if "_" in param:
                pname = param.split("_")[0]
                mname = "_".join(param.split("_")[1:])
                if sec in pbysec:
                    if mname not in pbysec[sec]:
                        pbysec[sec][mname] = {pname: v}
                    else:
                        pbysec[sec][mname][pname] = v
                else:
                    pbysec[sec] = {mname: {pname: v}}
            else:
                if sec in geom:
                    geom[sec][param] = v
                else:
                    geom[sec] = {param: v}
        return pbysec, geom

    PT5Bcells = {}
    if cfg.control:
        label = "control"
        plist = glob.glob("./cells/PT5B/Control*")
        vinit = vinit["Control animal PT neuron"]
        for fname in plist:
            if fname[13:-7] not in cfg.outliers:
                name = fname[38:-7]
                PT5Bcells[name] = {"file": fname}
    else:
        label = "6OHDA"
        plist = glob.glob("./cells/PT5B/6-OHDA*")
        vinit = vinit["6-OHDA animal PT neuron"]
        for fname in plist:
            if fname[13:-7] not in cfg.outliers:
                name = fname[37:-7]
                PT5Bcells[name] = {"file": fname}

    for i, name in enumerate(PT5Bcells):
        filename = PT5Bcells[name]["file"]
        # skip cells which did not scale well to the new morphology
        params, geom = paramsBySec(json.load(open(filename, "r")))
        fname = (
            f"Control animal PT neuron_{name}"
            if cfg.control
            else f"6-OHDA animal PT neuron_{name}"
        )
        rmp = vinit[name]
        cell = copy.deepcopy(PT5B.todict())
        for secname in cell["secs"]:
            cell["secs"][secname]["vinit"] = rmp
            sl = secname.split("_")[0] if "_" in secname else secname
            if "all" in geom:
                for k, v in geom["all"].items():
                    cell["secs"][secname]["geom"][k] = v
            if sl in geom:
                for k, v in geom[sl].items():
                    cell["secs"][secname]["geom"][k] = v
            if "all" in params:
                for k, v in params["all"].items():
                    cell["secs"][secname]["mechs"][k] = copy.copy(v)
            if sl in params:
                for k, v in params[sl].items():
                    cell["secs"][secname]["mechs"][k] = v
            cell["secs"][secname]["geom"]["nseg"] = cfg.pt5bNseg[secname]
            if cfg.removeNa:
                for naMech in cfg.naMechs:
                    if naMech in cell["secs"][secname]["mechs"]:
                        del cell["secs"][secname]["mechs"][naMech]
        cell["conds"]["cellType"] = f"PT5B{i}"
        # cell["conds"]["modelFit"] = name
        netParams.cellParams[f"PT5B{i}"] = cell
        if cfg.weightNorm:
            netParams.addCellParamsWeightNorm(
                f"PT5B{i}",
                f"conn/PT5B/PT5B{i}_{label}_weightNorm.pkl",
                threshold=cfg.weightNormThreshold,
            )
    del netParams.cellParams["PT5B"]


# ------------------------------------------------------------------------------
## IT5A full cell model params (700+ comps)
if "IT5A_full" not in loadCellParams:
    cellRule = netParams.importCellParams(
        label="IT5A_full",
        conds={"cellType": "IT", "cellModel": "HH_full", "ynorm": layer["5A"]},
        fileName="cells/ITcell.py",
        cellName="ITcell",
        cellArgs={"params": "BS1579"},
        somaAtOrigin=True,
    )
    netParams.renameCellParamsSec(label="IT5A_full", oldSec="soma_0", newSec="soma")
    netParams.addCellParamsWeightNorm(
        "IT5A_full",
        "conn/IT_full_BS1579_weightNorm.pkl",
        threshold=cfg.weightNormThreshold,
    )  # add weightNorm before renaming soma_0
    netParams.addCellParamsSecList(
        label="IT5A_full", secListName="perisom", somaDist=[0, 50]
    )  # sections within 50 um of soma
    cellRule["secLists"]["alldend"] = [
        sec for sec in cellRule.secs if ("dend" in sec or "apic" in sec)
    ]  # basal+apical
    cellRule["secLists"]["apicdend"] = [
        sec for sec in cellRule.secs if ("apic" in sec)
    ]  # basal+apical
    cellRule["secLists"]["spiny"] = [
        sec
        for sec in cellRule["secLists"]["alldend"]
        if sec not in ["apic_0", "apic_1"]
    ]
    # Adapt K gbar
    for secName in cellRule["secs"]:
        for kmech in [
            k
            for k in cellRule["secs"][secName]["mechs"].keys()
            if k.startswith("k") and k != "kBK"
        ]:
            cellRule["secs"][secName]["mechs"][kmech]["gbar"] *= cfg.KgbarFactor
    if saveCellParams:
        netParams.saveCellParamsRule(
            label="IT5A_full", fileName="cells/IT5A_full_cellParams.pkl"
        )


# ------------------------------------------------------------------------------
## IT5B full cell model params (700+ comps) - not used
# if 'IT5B_full' not in loadCellParams:
#   cellRule = netParams.importCellParams(label='IT5B_full', conds={'cellType': 'IT', 'cellModel': 'HH_full', 'ynorm': layer['5B']},
#     fileName='cells/ITcell.py', cellName='ITcell', cellArgs={'params': 'BS1579'}, somaAtOrigin=True)
#   netParams.addCellParamsSecList(label='IT5B_full', secListName='perisom', somaDist=[0, 50])  # sections within 50 um of soma
#   cellRule['secLists']['alldend'] = [sec for sec in cellRule.secs if ('dend' in sec or 'apic' in sec)] # basal+apical
#   cellRule['secLists']['apicdend'] = [sec for sec in cellRule.secs if ('apic' in sec)] # basal+apical
#   cellRule['secLists']['spiny'] = [sec for sec in cellRule['secLists']['alldend'] if sec not in ['apic_0', 'apic_1']]
#   netParams.addCellParamsWeightNorm('IT5B_full', 'conn/IT_full_BS1579_weightNorm.pkl')
#   netParams.saveCellParamsRule(label='IT5B_full', fileName='cells/IT5B_full_cellParams.pkl')


# ------------------------------------------------------------------------------
## PV cell params (3-comp)
if "PV_simple" not in loadCellParams:
    cellRule = netParams.importCellParams(
        label="PV_simple",
        conds={"cellType": "PV", "cellModel": "HH_simple"},
        fileName="cells/FS3.hoc",
        cellName="FScell1",
        cellInstance=True,
    )
    cellRule["secLists"]["spiny"] = ["soma", "dend"]
    netParams.addCellParamsWeightNorm(
        "PV_simple", "conn/PV_simple_weightNorm.pkl", threshold=cfg.weightNormThreshold
    )
    if saveCellParams:
        netParams.saveCellParamsRule(
            label="PV_simple", fileName="cells/PV_simple_cellParams.pkl"
        )


# ------------------------------------------------------------------------------
## SOM cell params (3-comp)
if "SOM_simple" not in loadCellParams:
    cellRule = netParams.importCellParams(
        label="SOM_simple",
        conds={"cellType": "SOM", "cellModel": "HH_simple"},
        fileName="cells/LTS3.hoc",
        cellName="LTScell1",
        cellInstance=True,
    )
    cellRule["secLists"]["spiny"] = ["soma", "dend"]
    netParams.addCellParamsWeightNorm(
        "SOM_simple",
        "conn/SOM_simple_weightNorm.pkl",
        threshold=cfg.weightNormThreshold,
    )
    if saveCellParams:
        netParams.saveCellParamsRule(
            label="SOM_simple", fileName="cells/SOM_simple_cellParams.pkl"
        )


# ------------------------------------------------------------------------------
# Population parameters
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
## load densities
with open("cells/cellDensity.pkl", "rb") as fileObj:
    density = pickle.load(fileObj)["density"]

## Local populations

netParams.popParams["IT2"] = {
    "cellModel": cfg.cellmod["IT2"],
    "cellType": "IT",
    "ynormRange": layer["2"],
    "density": density[("M1", "E")][0],
}
netParams.popParams["SOM2"] = {
    "cellModel": "HH_simple",
    "cellType": "SOM",
    "ynormRange": layer["24"],
    "density": density[("M1", "SOM")][5],
}
netParams.popParams["PV2"] = {
    "cellModel": "HH_simple",
    "cellType": "PV",
    "ynormRange": layer["24"],
    "density": density[("M1", "PV")][5],
}
netParams.popParams["IT4"] = {
    "cellModel": cfg.cellmod["IT4"],
    "cellType": "IT",
    "ynormRange": layer["4"],
    "density": density[("M1", "E")][1],
}
netParams.popParams["IT5A"] = {
    "cellModel": cfg.cellmod["IT5A"],
    "cellType": "IT",
    "ynormRange": layer["5A"],
    "density": density[("M1", "E")][2],
}
netParams.popParams["SOM5A"] = {
    "cellModel": "HH_simple",
    "cellType": "SOM",
    "ynormRange": layer["5A"],
    "density": density[("M1", "SOM")][2],
}
netParams.popParams["PV5A"] = {
    "cellModel": "HH_simple",
    "cellType": "PV",
    "ynormRange": layer["5A"],
    "density": density[("M1", "PV")][2],
}
netParams.popParams["IT5B"] = {
    "cellModel": cfg.cellmod["IT5B"],
    "cellType": "IT",
    "ynormRange": layer["5B"],
    "density": 0.5 * density[("M1", "E")][3],
}
if cfg.control is None:
    netParams.popParams["PT5B"] = {
        "cellModel": cfg.cellmod["PT5B"],
        "cellType": "PT",
        "ynormRange": layer["5B"],
        "density": 0.5 * density[("M1", "E")][3],
    }
else:
    for i, name in enumerate(PT5Bcells):
        netParams.popParams[f"PT5B{i}"] = {
            "cellType": f"PT5B{i}",
            "cellModel": "HH_reduced",
            "ynormRange": layer["5B"],
            "density": 0.5 * density[("M1", "E")][3] / len(PT5Bcells),
        }

netParams.popParams["SOM5B"] = {
    "cellModel": "HH_simple",
    "cellType": "SOM",
    "ynormRange": layer["5B"],
    "density": density[("M1", "SOM")][3],
}
netParams.popParams["PV5B"] = {
    "cellModel": "HH_simple",
    "cellType": "PV",
    "ynormRange": layer["5B"],
    "density": density[("M1", "PV")][3],
}
netParams.popParams["IT6"] = {
    "cellModel": cfg.cellmod["IT6"],
    "cellType": "IT",
    "ynormRange": layer["6"],
    "density": 0.5 * density[("M1", "E")][4],
}
netParams.popParams["CT6"] = {
    "cellModel": cfg.cellmod["CT6"],
    "cellType": "CT",
    "ynormRange": layer["6"],
    "density": 0.5 * density[("M1", "E")][4],
}
netParams.popParams["SOM6"] = {
    "cellModel": "HH_simple",
    "cellType": "SOM",
    "ynormRange": layer["6"],
    "density": density[("M1", "SOM")][4],
}
netParams.popParams["PV6"] = {
    "cellModel": "HH_simple",
    "cellType": "PV",
    "ynormRange": layer["6"],
    "density": density[("M1", "PV")][4],
}

if cfg.singleCellPops:
    for pop in netParams.popParams.values():
        pop["numCells"] = 1

# ------------------------------------------------------------------------------
## Long-range input populations (VecStims)
if cfg.addLongConn:
    ## load experimentally based parameters for long range inputs
    with open("conn/conn_long.pkl", "rb") as fileObj:
        connLongData = pickle.load(fileObj)
    # ratesLong = connLongData['rates']

    numCells = cfg.numCellsLong
    noise = cfg.noiseLong
    start = cfg.startLong

    longPops = ["TPO", "TVL", "S1", "S2", "cM1", "M2", "OC"]
    ## create populations with fixed
    for longPop in longPops:
        netParams.popParams[longPop] = {
            "cellModel": "VecStim",
            "numCells": numCells,
            "rate": cfg.ratesLong[longPop],
            "noise": noise,
            "start": start,
            "pulses": [],
            "ynormRange": layer["long" + longPop],
        }
        if isinstance(cfg.ratesLong[longPop], str):  # filename to load spikes from
            spikesFile = cfg.ratesLong[longPop]
            with open(spikesFile, "r") as f:
                spks = json.load(f)
            netParams.popParams[longPop].pop("rate")
            netParams.popParams[longPop]["spkTimes"] = spks


# ------------------------------------------------------------------------------
# Synaptic mechanism parameters
# ------------------------------------------------------------------------------
netParams.synMechParams["NMDA"] = {
    "mod": "MyExp2SynNMDABB",
    "tau1NMDA": 15,
    "tau2NMDA": 150,
    "e": 0,
}
netParams.synMechParams["AMPA"] = {
    "mod": "MyExp2SynBB",
    "tau1": 0.05,
    "tau2": 5.3 * cfg.AMPATau2Factor,
    "e": 0,
}
netParams.synMechParams["GABAB"] = {
    "mod": "MyExp2SynBB",
    "tau1": 3.5,
    "tau2": 260.9,
    "e": -93,
}
netParams.synMechParams["GABAA"] = {
    "mod": "MyExp2SynBB",
    "tau1": 0.07,
    "tau2": 18.2,
    "e": -80,
}
netParams.synMechParams["GABAASlow"] = {
    "mod": "MyExp2SynBB",
    "tau1": 2,
    "tau2": 100,
    "e": -80,
}
netParams.synMechParams["GABAASlowSlow"] = {
    "mod": "MyExp2SynBB",
    "tau1": 200,
    "tau2": 400,
    "e": -80,
}

ESynMech = ["AMPA", "NMDA"]
SOMESynMech = ["GABAASlow", "GABAB"]
SOMISynMech = ["GABAASlow"]
PVSynMech = ["GABAA"]


# ------------------------------------------------------------------------------
# Long range input pulses
# ------------------------------------------------------------------------------
if cfg.addPulses:
    for key in [k for k in dir(cfg) if k.startswith("pulse")]:
        params = getattr(cfg, key, None)
        [pop, start, end, rate, noise] = [
            params[s] for s in ["pop", "start", "end", "rate", "noise"]
        ]
        if (
            "duration" in params
            and params["duration"] is not None
            and params["duration"] > 0
        ):
            end = start + params["duration"]

        if pop in netParams.popParams:
            if "pulses" not in netParams.popParams[pop]:
                netParams.popParams[pop]["pulses"] = {}
            netParams.popParams[pop]["pulses"].append(
                {"start": start, "end": end, "rate": rate, "noise": noise}
            )


# ------------------------------------------------------------------------------
# Current inputs (IClamp)
# ------------------------------------------------------------------------------
if cfg.addIClamp:
    for key in [k for k in dir(cfg) if k.startswith("IClamp")]:
        params = getattr(cfg, key, None)
        [pop, sec, loc, start, dur, amp] = [
            params[s] for s in ["pop", "sec", "loc", "start", "dur", "amp"]
        ]

        # cfg.analysis['plotTraces']['include'].append((pop,0))  # record that pop

        # add stim source
        netParams.stimSourceParams[key] = {
            "type": "IClamp",
            "delay": start,
            "dur": dur,
            "amp": amp,
        }

        # connect stim source to target
        netParams.stimTargetParams[key + "_" + pop] = {
            "source": key,
            "conds": {"pop": pop},
            "sec": sec,
            "loc": loc,
        }

# ------------------------------------------------------------------------------
# NetStim inputs
# ------------------------------------------------------------------------------
if cfg.addNetStim:
    for key in [k for k in dir(cfg) if k.startswith("NetStim")]:
        params = getattr(cfg, key, None)
        [
            pop,
            ynorm,
            sec,
            loc,
            synMech,
            synMechWeightFactor,
            start,
            interval,
            noise,
            number,
            weight,
            delay,
        ] = [
            params[s]
            for s in [
                "pop",
                "ynorm",
                "sec",
                "loc",
                "synMech",
                "synMechWeightFactor",
                "start",
                "interval",
                "noise",
                "number",
                "weight",
                "delay",
            ]
        ]

        # cfg.analysis['plotTraces']['include'] = [(pop,0)]

        if synMech == ESynMech:
            wfrac = cfg.synWeightFractionEE
        elif synMech == SOMESynMech:
            wfrac = cfg.synWeightFractionSOME
        else:
            wfrac = [1.0]

        # add stim source
        netParams.stimSourceParams[key] = {
            "type": "NetStim",
            "start": start,
            "interval": interval,
            "noise": noise,
            "number": number,
        }

        # connect stim source to target
        # for i, syn in enumerate(synMech):
        netParams.stimTargetParams[key + "_" + pop] = {
            "source": key,
            "conds": {"pop": pop, "ynorm": ynorm},
            "sec": sec,
            "loc": loc,
            "synMech": synMech,
            "weight": weight,
            "synMechWeightFactor": synMechWeightFactor,
            "delay": delay,
        }

# ------------------------------------------------------------------------------
# Local connectivity parameters
# ------------------------------------------------------------------------------
with open("conn/conn.pkl", "rb") as fileObj:
    connData = pickle.load(fileObj)
pmat = connData["pmat"]
wmat = connData["wmat"]
bins = connData["bins"]
# import conn_param
# pmat = conn_param.pmat
# wmat = conn_param.wmat
# bins = conn_param.bins

# ------------------------------------------------------------------------------
## E -> E
if cfg.addConn and cfg.EEGain > 0.0:
    PTs = ["PT"] if cfg.control is None else [f"PT5B{i}" for i in range(len(PT5Bcells))]
    labelsConns = [
        ("W+AS_norm", "IT", "L2/3,4"),
        ("W+AS_norm", "IT", "L5A,5B"),
        ("W+AS_norm", "PT", "L5B"),
        ("W+AS_norm", "IT", "L6"),
        ("W+AS_norm", "CT", "L6"),
    ]
    labelPostBins = [
        ("W+AS", "IT", "L2/3,4"),
        ("W+AS", "IT", "L5A,5B"),
        ("W+AS", "PT", "L5B"),
        ("W+AS", "IT", "L6"),
        ("W+AS", "CT", "L6"),
    ]
    labelPreBins = ["W", "AS", "AS", "W", "W"]
    preTypes = [["IT"], ["IT"], ["IT"] + PTs, ["IT", "CT"], ["IT", "CT"]]
    postTypes = [["IT"], ["IT"], PTs, ["IT"], ["CT"]]
    ESynMech = ["AMPA", "NMDA"]

    for i, (label, preBinLabel, postBinLabel) in enumerate(
        zip(labelsConns, labelPreBins, labelPostBins)
    ):
        for ipre, preBin in enumerate(bins[preBinLabel]):
            for ipost, postBin in enumerate(bins[postBinLabel]):
                for cellModel in cellModels:
                    ruleLabel = (
                        "EE_"
                        + cellModel
                        + "_"
                        + str(i)
                        + "_"
                        + str(ipre)
                        + "_"
                        + str(ipost)
                    )
                    netParams.connParams[ruleLabel] = {
                        "preConds": {"cellType": preTypes[i], "ynorm": list(preBin)},
                        "postConds": {
                            "cellModel": cellModel,
                            "cellType": postTypes[i],
                            "ynorm": list(postBin),
                        },
                        "synMech": ESynMech,
                        "probability": pmat[label][ipost, ipre],
                        "weight": wmat[label][ipost, ipre]
                        * cfg.EEGain
                        / cfg.synsperconn[cellModel],
                        "synMechWeightFactor": cfg.synWeightFractionEE,
                        "delay": "defaultDelay+dist_3D/propVelocity",
                        "synsPerConn": cfg.synsperconn[cellModel],
                        "sec": "spiny",
                    }


# ------------------------------------------------------------------------------
## E -> I
if cfg.EIGain:  # Use IEGain if value set
    cfg.EPVGain = cfg.EIGain
    cfg.ESOMGain = cfg.EIGain
else:
    cfg.EIGain = (cfg.EPVGain + cfg.ESOMGain) / 2.0

if cfg.addConn and (cfg.EPVGain > 0.0 or cfg.ESOMGain > 0.0):
    labelsConns = ["FS", "LTS"]
    labelPostBins = ["FS/LTS", "FS/LTS"]
    labelPreBins = ["FS/LTS", "FS/LTS"]
    preTypes = ["IT", "PT", "CT"]
    postTypes = ["PV", "SOM"]
    ESynMech = ["AMPA", "NMDA"]
    lGain = [cfg.EPVGain, cfg.ESOMGain]  # E -> PV or E -> SOM
    for i, (label, preBinLabel, postBinLabel) in enumerate(
        zip(labelsConns, labelPreBins, labelPostBins)
    ):
        for ipre, preBin in enumerate(bins[preBinLabel]):
            for ipost, postBin in enumerate(bins[postBinLabel]):
                ruleLabel = "EI_" + str(i) + "_" + str(ipre) + "_" + str(ipost)
                netParams.connParams[ruleLabel] = {
                    "preConds": {"cellType": preTypes, "ynorm": list(preBin)},
                    "postConds": {"cellType": postTypes[i], "ynorm": list(postBin)},
                    "synMech": ESynMech,
                    "probability": pmat[label][ipost, ipre],
                    "weight": wmat[label][ipost, ipre] * lGain[i],
                    "synMechWeightFactor": cfg.synWeightFractionEI,
                    "delay": "defaultDelay+dist_3D/propVelocity",
                    "sec": "soma",
                }  # simple I cells used right now only have soma


# ------------------------------------------------------------------------------
## I -> all
if cfg.IEGain:  # Use IEGain if value set
    cfg.PVEGain = cfg.IEGain
    cfg.SOMEGain = cfg.IEGain
else:
    cfg.IEGain = (cfg.PVEGain + cfg.SOMEGain) / 2.0

if cfg.IIGain:  # Use IIGain if value set
    cfg.SOMPVGain = cfg.IIGain
    cfg.PVSOMGain = cfg.IIGain
    cfg.SOMSOMGain = cfg.IIGain
    cfg.PVPVGain = cfg.IIGain
else:
    cfg.IIGain = (cfg.PVSOMGain + cfg.SOMPVGain + cfg.SOMSOMGain + cfg.PVPVGain) / 4.0

if cfg.addConn and (cfg.IEGain > 0.0 or cfg.IIGain > 0.0):
    # Local, intralaminar only; all-to-all but distance-based; high weights; L5A/B->L5A/B
    preCellTypes = ["SOM", "SOM", "SOM", "PV", "PV", "PV"]
    ynorms = [
        (layer["2"][0], layer["4"][1]),
        (layer["5A"][0], layer["5B"][1]),
        (layer["6"][0], layer["6"][1]),
    ] * 2
    IEweights = (
        cfg.IEweights * 2
    )  # [I->E2/3+4, I->E5, I->E6] weights (Note * 2 is repeat list operator)
    IIweights = (
        cfg.IIweights * 2
    )  # [I->I2/3+4, I->I5, I->I6] weights (Note * 2 is repeat list operator)
    postCellTypes = ["PT", ["IT", "CT"], "PV", "SOM"]
    IEdisynBiases = [
        None,
        cfg.IEdisynapticBias,
        cfg.IEdisynapticBias,
        None,
        cfg.IEdisynapticBias,
        cfg.IEdisynapticBias,
    ]
    disynapticBias = None  # default, used for I->I

    # preCellTypes = ['SOM', 'PV']
    # ynorms = [[0,1]]*2
    # IEweights = cfg.IEweights * 2  # [I->E2/3+4, I->E5, I->E6] weights (Note * 2 is repeat list operator)
    # IIweights = cfg.IIweights * 2  # [I->I2/3+4, I->I5, I->I6] weights (Note * 2 is repeat list operator)
    # postCellTypes = ['PT', ['IT','CT'], 'PV', 'SOM']
    # IEdisynBiases = [cfg.IEdisynapticBias, cfg.IEdisynapticBias]
    # disynapticBias = None  # default, used for I->I

    for i, (preCellType, ynorm, IEweight, IIweight, IEdisynBias) in enumerate(
        zip(preCellTypes, ynorms, IEweights, IIweights, IEdisynBiases)
    ):
        for ipost, postCellType in enumerate(postCellTypes):
            for cellModel in cellModels:
                if postCellType == "PV":  # postsynaptic I cell
                    sec = "soma"
                    synWeightFraction = [1]
                    if preCellType == "PV":  # PV->PV
                        weight = IIweight * cfg.PVPVGain
                        synMech = PVSynMech
                    else:  # SOM->PV
                        weight = IIweight * cfg.SOMPVGain
                        synMech = SOMISynMech
                elif postCellType == "SOM":  # postsynaptic I cell
                    sec = "soma"
                    synWeightFraction = [1]
                    if preCellType == "PV":  # PV->SOM
                        weight = IIweight * cfg.PVSOMGain
                        synMech = PVSynMech
                    else:  # SOM->SOM
                        weight = IIweight * cfg.SOMSOMGain
                        synMech = SOMISynMech
                elif postCellType == ["IT", "CT"]:  # postsynaptic IT,CT cell
                    disynapticBias = IEdisynBias
                    if preCellType == "PV":  # PV->E
                        weight = IEweight * cfg.PVEGain
                        synMech = PVSynMech
                        sec = "perisom"
                    else:  # SOM->E
                        weight = IEweight * cfg.SOMEGain
                        synMech = SOMESynMech
                        sec = "spiny"
                        synWeightFraction = cfg.synWeightFractionSOME
                elif postCellType == "PT":  # postsynaptic PT cell
                    disynapticBias = IEdisynBias
                    if preCellType == "PV":  # PV->E
                        weight = IEweight * cfg.IPTGain * cfg.PVEGain
                        synMech = PVSynMech
                        sec = "perisom"
                    else:  # SOM->E
                        weight = IEweight * cfg.IPTGain * cfg.SOMEGain
                        synMech = SOMESynMech
                        sec = "spiny"
                        synWeightFraction = cfg.synWeightFractionSOME
                if cellModel == "HH_full":
                    weight = weight * cfg.IFullGain

                ruleLabel = "I_" + cellModel + "_" + str(i) + "_" + str(ipost)
                netParams.connParams[ruleLabel] = {
                    "preConds": {"cellType": preCellType, "ynorm": ynorm},
                    "postConds": {
                        "cellModel": cellModel,
                        "cellType": postCellType,
                        "ynorm": ynorm,
                    },
                    "synMech": synMech,
                    "probability": "1.0 * exp(-dist_3D_border/probLambda)",
                    "weight": weight / cfg.synsperconn[cellModel],
                    "delay": "defaultDelay+dist_3D_border/propVelocity",
                    "synsPerConn": cfg.synsperconn[cellModel],
                    "synMechWeightFactor": synWeightFraction,
                    "sec": sec,
                    "disynapticBias": disynapticBias,
                }


# ------------------------------------------------------------------------------
# Long-range  connectivity parameters
# ------------------------------------------------------------------------------
if cfg.addLongConn:
    # load load experimentally based parameters for long range inputs
    cmatLong = connLongData["cmat"]
    binsLong = connLongData["bins"]

    longPops = ["TPO", "TVL", "S1", "S2", "cM1", "M2", "OC"]
    PTs = ["PT"] if cfg.control is None else [f"PT5B{i}" for i in range(len(PT5Bcells))]
    cellTypes = ["IT"] + PTs + ["CT", "PV", "SOM"]
    cellTypeLabels = ["IT"] + ["PT" for _ in PTs] + ["CT", "PV", "SOM"]
    EorI = ["exc", "inh"]
    syns = {"exc": ESynMech, "inh": "GABAA"}
    synFracs = {"exc": cfg.synWeightFractionEE, "inh": [1.0]}

    for longPop in longPops:
        for ct, lab in zip(cellTypes, cellTypeLabels):
            for EorI in ["exc", "inh"]:
                for i, (binRange, convergence) in enumerate(
                    zip(binsLong[(longPop, lab)], cmatLong[(longPop, lab, EorI)])
                ):
                    for cellModel in cellModels:
                        ruleLabel = (
                            longPop
                            + "_"
                            + ct
                            + "_"
                            + EorI
                            + "_"
                            + cellModel
                            + "_"
                            + str(i)
                        )
                        netParams.connParams[ruleLabel] = {
                            "preConds": {"pop": longPop},
                            "postConds": {
                                "cellModel": cellModel,
                                "cellType": ct,
                                "ynorm": list(binRange),
                            },
                            "synMech": syns[EorI],
                            "convergence": convergence,
                            "weight": cfg.weightLong / cfg.synsperconn[cellModel],
                            "synMechWeightFactor": cfg.synWeightFractionEE,
                            "delay": "defaultDelay+dist_3D/propVelocity",
                            "synsPerConn": (
                                cfg.synsperconn[lab]
                                if lab == "PT" and cfg.control is not None
                                else cfg.synsperconn[cellModel]
                            ),
                            "sec": "spiny",
                        }


# ------------------------------------------------------------------------------
# Subcellular connectivity (synaptic distributions)
# ------------------------------------------------------------------------------
PT5B = ["PT5B"] if cfg.control is None else [f"PT5B{i}" for i in range(len(PT5Bcells))]

if cfg.addSubConn:
    with open("conn/conn_dend_PT.json", "r") as fileObj:
        connDendPTData = json.load(fileObj)
    with open("conn/conn_dend_IT.json", "r") as fileObj:
        connDendITData = json.load(fileObj)

    # ------------------------------------------------------------------------------
    # L2/3,TVL,S2,cM1,M2 -> PT (Suter, 2015)
    lenY = 30
    spacing = 50
    gridY = list(range(0, -spacing * lenY, -spacing))
    synDens, _, fixedSomaY = (
        connDendPTData["synDens"],
        connDendPTData["gridY"],
        connDendPTData["fixedSomaY"],
    )
    for k in synDens.keys():
        prePop, postLabel = k.split("_")  # eg. split 'M2_PT'
        if prePop == "L2":
            prePop = "IT2"  # include conns from layer 2/3 and 4

        # iterate over different PT5B cell types
        if postLabel == "PT":
            postTypes = PT5B
        else:
            postTypes = [postLabel]
        for postType in postTypes:
            netParams.subConnParams[k] = {
                "preConds": {"pop": prePop},
                "postConds": {"cellType": postType},
                "sec": "spiny",
                "groupSynMechs": ESynMech,
                "density": {
                    "type": "1Dmap",
                    "gridX": None,
                    "gridY": gridY,
                    "gridValues": synDens[k],
                    "fixedSomaY": fixedSomaY,
                },
            }

    # ------------------------------------------------------------------------------
    # TPO, TVL, M2, OC  -> E (L2/3, L5A, L5B, L6) (Hooks 2013)
    lenY = 26
    spacing = 50
    gridY = list(range(0, -spacing * lenY, -spacing))
    synDens, _, fixedSomaY = (
        connDendITData["synDens"],
        connDendITData["gridY"],
        connDendITData["fixedSomaY"],
    )
    for k in synDens.keys():
        prePop, post = k.split("_")  # eg. split 'M2_L2'
        postCellTypes = (
            ["IT", "PT", "CT"] if prePop in ["OC", "TPO"] else ["IT", "CT"]
        )  # only OC,TPO include PT cells
        postyRange = list(layer[post.split("L")[1]])  # get layer yfrac range
        if post == "L2":
            postyRange[1] = layer["4"][1]  # apply L2 rule also to L4
        netParams.subConnParams[k] = {
            "preConds": {"pop": prePop},
            "postConds": {"ynorm": postyRange, "cellType": postCellTypes},
            "sec": "spiny",
            "groupSynMechs": ESynMech,
            "density": {
                "type": "1Dmap",
                "gridX": None,
                "gridY": gridY,
                "gridValues": synDens[k],
                "fixedSomaY": fixedSomaY,
            },
        }

    # ------------------------------------------------------------------------------
    # S1, S2, cM1 -> E IT/CT; no data, assume uniform over spiny
    netParams.subConnParams["S1,S2,cM1->IT,CT"] = {
        "preConds": {"pop": ["S1", "S2", "cM1"]},
        "postConds": {"cellType": ["IT", "CT"]},
        "sec": "spiny",
        "groupSynMechs": ESynMech,
        "density": "uniform",
    }

    # ------------------------------------------------------------------------------
    # rest of local E->E (exclude IT2->PT); uniform distribution over spiny
    netParams.subConnParams["IT2->non-PT"] = {
        "preConds": {"pop": ["IT2"]},
        "postConds": {"cellType": ["IT", "CT"]},
        "sec": "spiny",
        "groupSynMechs": ESynMech,
        "density": "uniform",
    }

    netParams.subConnParams["non-IT2->E"] = {
        "preConds": {"pop": ["IT4", "IT5A", "IT5B", "IT6", "CT6"] + PT5B},
        "postConds": {"cellType": ["IT", "PT", "CT"]},
        "sec": "spiny",
        "groupSynMechs": ESynMech,
        "density": "uniform",
    }

    # ------------------------------------------------------------------------------
    # PV->E; perisomatic (no sCRACM)
    netParams.subConnParams["PV->E"] = {
        "preConds": {"cellType": "PV"},
        "postConds": {"cellType": ["IT", "CT", "PT"]},
        "sec": "perisom",
        "density": "uniform",
    }

    # ------------------------------------------------------------------------------
    # SOM->E; apical dendrites (no sCRACM)
    netParams.subConnParams["SOM->E"] = {
        "preConds": {"cellType": "SOM"},
        "postConds": {"cellType": ["IT", "CT", "PT"]},
        "sec": "apicdend",
        "groupSynMechs": SOMESynMech,
        "density": "uniform",
    }

    # ------------------------------------------------------------------------------
    # All->I; apical dendrites (no sCRACM)
    netParams.subConnParams["All->I"] = {
        "preConds": {"cellType": ["IT", "CT", "PT", "SOM", "PV"] + longPops},
        "postConds": {"cellType": ["SOM", "PV"]},
        "sec": "spiny",
        "groupSynMechs": ESynMech,
        "density": "uniform",
    }


# ------------------------------------------------------------------------------
# Description
# ------------------------------------------------------------------------------
netParams.description = """ 
- M1 net, 6 layers, 7 cell types 
- NCD-based connectivity from  Weiler et al. 2008; Anderson et al. 2010; Kiritani et al. 2012; 
  Yamawaki & Shepherd 2015; Apicella et al. 2012
- Parametrized version based on Sam's code
- Updated cell models and mod files
- Added parametrized current inputs
- Fixed bug: prev was using cell models in /usr/site/nrniv/local/python/ instead of cells 
- Use 5 synsperconn for 5-comp cells (HH_reduced); and 1 for 1-comp cells (HH_simple)
- Fixed bug: made global h params separate for each cell model
- Fixed v_init for different cell models
- New IT cell with same geom as PT
- Cleaned cfg and moved background inputs here
- Set EIGain and IEGain for each inh cell type
- Added secLists for PT full
- Fixed reduced CT (wrong vinit and file)
- Added subcellular conn rules to distribute synapses
- PT full model soma centered at 0,0,0 
- Set cfg seeds here to ensure they get updated
- Added PVSOMGain and SOMPVGain
- PT subcellular distribution as a cfg param
- Cylindrical volume
- DefaultDelay (for local conns) = 2ms
- Added long range connections based on Yamawaki 2015a,b; Suter 2015; Hooks 2013; Meyer 2011
- Updated cell densities based on Tsai 2009; Lefort 2009; Katz 2011; Wall 2016; 
- Separated PV and SOM of L5A vs L5B
- Fixed bugs in local conn (PT, PV5, SOM5, L6)
- Added perisom secList including all sections 50um from soma
- Added subcellular conn rules (for both full and reduced models)
- Improved cell models, including PV and SOM fI curves
- Improved subcell conn rules based on data from Suter15, Hooks13 and others
- Adapted Bdend L of reduced cell models
- Made long pop rates a cfg param
- Set threshold to 0.0 mV
- Parametrized I->E/I layer weights
- Added missing subconn rules (IT6->PT; S1,S2,cM1->IT/CT; long->SOM/PV)
- Added threshold to weightNorm (PT threshold=10x)
- weightNorm threshold as a cfg parameter
- Separate PV->SOM, SOM->PV, SOM->SOM, PV->PV gains 
- Conn changes: reduced IT2->IT4, IT5B->CT6, IT5B,6->IT2,4,5A, IT2,4,5A,6->IT5B; increased CT->PV6+SOM6
- Parametrized PT ih gbar
- Added IFullGain parameter: I->E gain for full detailed cell models
- Replace PT ih with Migliore 2012
- Parametrized ihGbar, ihGbarBasal, dendNa, axonNa, axonRa, removeNa
- Replaced cfg list params with dicts
- Parametrized ihLkcBasal and AMPATau2Factor
- Fixed synMechWeightFactor
- Parametrized PT ih slope
- Added disynapticBias to I->E (Yamawaki&Shepherd,2015)
- Fixed E->CT bin 0.9-1.0
- Replaced GABAB with exp2syn and adapted synMech ratios
- Parametrized somaNa
- Added ynorm condition to NetStims
- Added option to play back recorded spikes into long-range inputs
- Fixed Bdend pt3d y location
- Added netParams.convertCellShapes = True to convert stylized geoms to 3d points
- New layer boundaries, cell densities, conn, FS+SOM L4 grouped with L2/3, low cortical input to L4
- Increased exc->L4 based on Yamawaki 2015 fig 5
- v54: Moved from NetPyNE v0.7.9 to v0.9.1 (v54_batch1-6)
- v54: Moved to NetPyNE v0.9.1 and py3 (v54_batch7 onwards)
- v56: Reduced dt from 0.05 to 0.025 (note this version follows from v54, i.e. without new cell types; branch 'paper2019_py3')
- v56: (included in prev version): Added cfg.KgbarFactor
"""
