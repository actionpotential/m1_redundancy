"""
initWeightNorm.py

Find weights that give a 0.5 mV somatic response for the PT5B cells.

"""

import matplotlib

#matplotlib.use("Agg")  # to avoid graphics error in servers
from netpyne import sim
from neuron import h
import json
import pickle

# Load cfg.py and modify with single cells and no connections
from cfg import cfg

cfg.singleCellPops = True
cfg.weightNorm = False
cfg.addConn = False
cfg.addSubConn = False
cfg.addLongConn = False
cfg.addIClamp = False
cfg.addNetStim = True
cfg.duration = 200
cfg.dt = 0.01
cfg.recDt = cfg.dt
secs = ["soma", "dend_0", "apic_0"]

# record every section
cfg.recordTraces = {f"v_soma": {"sec": "soma", "loc": 0.5, "var": "v"}}

Ncell = cfg.NPTCells 
ns = {f"PT5B{pop}": {} for pop in range(Ncell)}
nsegs = [cfg.pt5bNseg[sec] for sec in secs]
bounds = {f"PT5B{pop}": {sec: [(0, 0.001) for _ in range(n)] for sec,n in zip(secs,nsegs)} for pop in range(Ncell)}

# Add a NetStim to each section of each PT5B neuron
i = 1
for sec,n in zip(secs,nsegs):
    for pop in range(Ncell):
        ns[f"PT5B{pop}"][sec] = []
        for j in range(n):
            ns[f"PT5B{pop}"][sec].append(f"NetStim{i}")
            lb, up = bounds[f"PT5B{pop}"][sec][j]
            cfg[f"NetStim{i}"] = {
                "pop": f"PT5B{pop}",
                "ynorm": [0, 1],
                "sec": "sec",
                "loc": (j+1)/(n+1),
                "synMech": ["AMPA"],
                "synMechWeightFactor": [1.0],
                "start": 150,
                "interval": 10,
                "noise": 0.0,
                "number": 1,
                "weight": 0,
                "delay": 0,
            }
            i += 1
start_idx = int(150 / 0.025)
# Load netParams which import cfg from __main__
from netParams import netParams
for pop in list(netParams.popParams.keys()):
    if pop not in bounds:
        del netParams.popParams[pop]

for sec,n in zip(secs,nsegs):
    for j in range(n):
        for p in range(Ncell):
            pop = f"PT5B{p}"
            lb, up = bounds[pop][sec][j]
            netParams.stimTargetParams[f"{ns[pop][sec][j]}_{pop}"]["weight"] = (lb + up) / 2
        sim.create(simConfig=cfg, netParams=netParams, clearAll=sec != secs[0] and j != 0)
        sim.runSim()
        for _ in range(30):
            for p in range(Ncell):
                pop = f"PT5B{p}"
                idx = sim.net.pops[pop].cellGids[0]
                v = sim.simData[f"v_soma"][f"cell_{idx}"].as_numpy()
                step = max(v[start_idx:]) - v[start_idx]
                lb, up = bounds[pop][sec][j] 
                weight = (lb + up) / 2
                if step > 0.5:
                    up = weight
                else:
                    lb = weight
                bounds[pop][sec][j] = (lb, up)
                netParams.stimTargetParams[f"{ns[pop][sec][j]}_{pop}"]["weight"] = (
                    lb + up
                ) / 2
                """
                modify = {'conds': {'preLabel': ns[pop][sec]},
                          'postConds': {'pop':pop},
                          'weight': (lb + up)/2}
                sim.net.modifyConns(modify)
                """
            # modifyConns slows down each simulation
            # TODO: find out why
            sim.create(simConfig=cfg, netParams=netParams, clearAll=True)
            sim.runSim()
        # reset weights
        for p in range(Ncell):
            pop = f"PT5B{p}"
            netParams.stimTargetParams[f"{ns[pop][sec][j]}_{pop}"]["weight"] = 0


# save the new weights
label = "control" if cfg.control else "6OHDA"
for pop in bounds:
    weights = {}
    for sec, bnd in bounds[pop].items():
        weights[sec] = [(u+l)/2 for (u,l) in bnd]
    pickle.dump(weights, open(f"conn/PT5B/{pop}_{label}_weightNorm.pkl", "wb"))
