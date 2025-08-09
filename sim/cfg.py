"""
cfg.py 

Simulation configuration for M1 model (using NetPyNE)

Contributors: salvadordura@gmail.com
"""

from netpyne import specs
from datetime import date
import pickle

cfg = specs.SimConfig()

# ------------------------------------------------------------------------------
#
# SIMULATION CONFIGURATION
#
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Run parameters
# ------------------------------------------------------------------------------
cfg.duration = 4300
cfg.dt = 0.025
cfg.seeds = {"conn": 2484, "stim": 9931, "loc": 8516}
cfg.hParams = {"celsius": 34, "v_init": -80}
cfg.verbose = False
cfg.createNEURONObj = 1
cfg.createPyStruct = 1
cfg.connRandomSecFromList = False  # set to false for reproducibility
cfg.cvode_active = False
cfg.cvode_atol = 1e-6
cfg.cache_efficient = False
cfg.printRunTime = 0.1
cfg.oneSynPerNetcon = (
    True  # only affects conns not in subconnParams; produces identical results
)

cfg.includeParamsLabel = False  # True # needed for modify synMech False
cfg.printPopAvgRates = [0, 400.0]

cfg.checkErrors = True

# select cell from healthy controls (True) or 6-OHDA mice (False)
cfg.control = False
cfg.pt5bEna = 53
cfg.pt5bEk = -104
cfg.pt5bNseg = {
    "soma": 1,
    "apic_0": 11,
    "dend_0": 11,
    "axon_0": 3,
}  # nseg for the apical and basal dendrites for the reduced model
# cfg.saveInterval = 100 # define how often the data is saved, this can be used with interval run if you want to update the weights more often than you save
# cfg.intervalFolder = 'interval_saving'


# ------------------------------------------------------------------------------
# Recording
# ------------------------------------------------------------------------------
cfg.NPTCells = 28 if cfg.control else 16
allpops = [
    "IT2",
    "PV2",
    "SOM2",
    "IT4",
    "IT5A",
    "PV5A",
    "SOM5A",
    "IT5B",
    "PV5B",
    "SOM5B",
    "IT6",
    "CT6",
    "PV6",
    "SOM6",
] + [f"PT5B{i}" for i in range(cfg.NPTCells)]
allPTpops = [f"PT5B{i}" for i in range(cfg.NPTCells)]

cfg.cellsrec = 1 
if cfg.cellsrec == 0:  cfg.recordCells = ['all'] # record all cells
elif cfg.cellsrec == 1: cfg.recordCells = [(pop, i) for pop in allPTpops for i in [9]] # record about 10 neurons from each population of PT neurons
elif cfg.cellsrec == 2: cfg.recordCells = [(pop, i) for pop in allpops for i in [0]]  # record from 1 cell of  of each population
elif cfg.cellsrec == 3: cfg.recordCells = [(pop,50) for pop in ['IT5A', 'PT5B']]+[('PT5B',x) for x in [393,579,19,104]] #,214,1138,799]] # record selected cells # record selected cells
elif cfg.cellsrec == 4: cfg.recordCells = [(pop,50) for pop in ['IT2', 'IT4', 'IT5A', 'PT5B']]+[('PT5B',x) for x in [393,447,579,19,104,214,1138,979,799]] # record selected cells

cfg.recordTraces = {"v_soma": {"sec": "soma", "loc": 0.5, "var": "v"}}  # ,
#'V_soma_ih': {'sec':'soma', 'loc':0.5, 'var':'gbar', 'mech':'hd', 'conds':{'pop': 'PT5B'}}}
# 'V_apic_26': {'sec':'apic_26', 'loc':0.5, 'var':'v', 'conds':{'pop': 'PT5B'}},
# 'V_dend_5': {'sec':'dend_5', 'loc':0.5, 'var':'v', 'conds':{'pop': 'PT5B'}}}
#'I_AMPA_Adend2': {'sec':'Adend2', 'loc':0.5, 'synMech': 'AMPA', 'var': 'i'}}

cfg.recordLFP = [
    [150, y, 150] for y in range(200, 1300, 100)
]  # [[150, y, 150] for y in range(200,1300,100)]

cfg.saveLFPPops = False  # allpops

cfg.recordDipoles = (
    False  # {'L2': ['IT2'], 'L4': ['IT4'], 'L5': ['IT5A', 'IT5B', 'PT5B']}
)

cfg.recordStim = False
cfg.recordTime = False
cfg.recordStep = 0.025
cfg.outliers = [
    "Control animal PT neuron_20201215-B6-2826_cell10",
    "Control animal PT neuron_20201215-B6-2826_cell3",
    "Control animal PT neuron_20201215-B6-2826_cell2",
    "Control animal PT neuron_20201230-B6-2894_cell7",
    "Control animal PT neuron_20201230-B6-2894_cell6",
    "Control animal PT neuron_20201230-B6-2894_cell5",
    "Control animal PT neuron_20201230-B6-2894_cell1",
    "Control animal PT neuron_20201201-03445-07258_cell3",
    "Control animal PT neuron_20201201-03445-07258_cell5",
    "Control animal PT neuron_20201210-B6-2825_cell9",
    "Control animal PT neuron_20201210-B6-2825_cell2",
] + [
    "6-OHDA animal PT neuron_20201202-03466-07377_cell1",
    "6-OHDA animal PT neuron_20201222-B6-2852_cell8",
    "6-OHDA animal PT neuron_20201222-B6-2852_cell6",
    "6-OHDA animal PT neuron_20201222-B6-2852_cell1",
    "6-OHDA animal PT neuron_20201222-B6-2852_cell7",
    "6-OHDA animal PT neuron_20201222-B6-2852_cell10",
    "6-OHDA animal PT neuron_20201222-B6-2852_cell3",
    "6-OHDA animal PT neuron_20201222-B6-2852_cell9",
    "6-OHDA animal PT neuron_20201211-03466-07376_cell7",
    "6-OHDA animal PT neuron_20201211-03466-07376_cell2",
    "6-OHDA animal PT neuron_20201211-03466-07376_cell3",
    "6-OHDA animal PT neuron_20201211-03466-07376_cell9",
]

# ------------------------------------------------------------------------------
# Saving
# ------------------------------------------------------------------------------
cfg.simLabel = "aM1_07-30-2025_01"
cfg.saveFolder = "data/"
cfg.savePickle = False
cfg.saveJson = True
cfg.saveDataInclude = ["simData", "simConfig", "netParams", "net"]
cfg.backupCfgFile = None  # ['cfg.py', 'backupcfg/']
cfg.gatherOnlySimData = False
cfg.saveCellSecs = False
cfg.saveCellConns = False
cfg.compactConnFormat = 0


# ------------------------------------------------------------------------------
# Analysis and plotting
# ------------------------------------------------------------------------------
with open("cells/popColors.pkl", "rb") as fileObj:
    popColors = pickle.load(fileObj)["popColors"]

#cfg.analysis["plotRaster"] = {
#    "include": allpops,
#    "orderBy": ["pop", "y"],
#    "timeRange": [0, cfg.duration],
#    "saveFig": True,
#    "showFig": False,
#    "popRates": True,
#    "orderInverse": True,
#    "popColors": popColors,
#    "figSize": (12, 10),
#    "lw": 0.3,
#    "markerSize": 3,
#    "marker": ".",
#    "dpi": 300,
#}


# cfg.analysis['plotLFP'] = {'plots': ['timeSeries'], 'electrodes': list(range(len(cfg.recordLFP))), 'figSize': (12,10), 'timeRange': [1000,5000],  'saveFig': True, 'showFig':False}


# cfg.analysis['plotTraces'] = {'include': [], 'timeRange': [0, cfg.duration], 'oneFigPer': 'trace', 'figSize': (10,4), 'saveFig': True, 'showFig': False}


# ------------------------------------------------------------------------------
# Cells
# ------------------------------------------------------------------------------
cfg.cellmod = {
    "IT2": "HH_reduced",
    "IT4": "HH_reduced",
    "IT5A": "HH_full",
    "IT5B": "HH_reduced",
    "PT5B": "HH_full" if cfg.control is None else "HH_reduced",
    "IT6": "HH_reduced",
    "CT6": "HH_reduced",
}

cfg.ihModel = "migliore"  # ih model
cfg.ihGbar = 0.75  # multiplicative factor for ih gbar in PT cells

cfg.ihGbarZD = None  # multiplicative factor for ih gbar in PT cells
cfg.ihGbarBasal = 1.0  # 0.1 # multiplicative factor for ih gbar in PT cells
cfg.ihlkc = 0.2  # ih leak param (used in Migliore)
cfg.ihlkcBasal = 1.0
cfg.ihlkcBelowSoma = 0.01
cfg.ihlke = -86  # ih leak param (used in Migliore)
cfg.ihSlope = 14 * 2

cfg.removeNa = False  # simulate TTX; set gnabar=0s

cfg.naMechs = ["Nap_Et2", "NaTs2_t", "NaTa_t", "Nap_Et2"]
cfg.somaNa = 5
cfg.dendNa = 0.3
cfg.axonNa = 7
cfg.axonRa = 0.005

cfg.gpas = 0.5  # multiplicative factor for pas g in PT cells
cfg.epas = 0.9  # multiplicative factor for pas e in PT cells

cfg.KgbarFactor = 1.0  # multiplicative factor for K channels gbar in all E cells
cfg.makeKgbarFactorEqualToNewFactor = False

cfg.modifyMechs = {
    "startTime": 500,
    "endTime": 1000,
    "cellType": "PT",
    "mech": "hd",
    "property": "gbar",
    "newFactor": 1.00,
    "origFactor": 0.75,
}

# ------------------------------------------------------------------------------
# Channel Parameters to be changed
# ------------------------------------------------------------------------------
cfg.gbar = 0.0055  # NaP from nap_sidi.mod; 0.0153130 Sam 2017 ... only works with rest standard when >= 0.0055
cfg.gpeak = 7.251280172010002e-05  # BK 7.251280172010002e-05 Sam 2017
cfg.gnafbar = 0.00086  # NaT from nafx.mod; 0.00086 control

# ------------------------------------------------------------------------------
# Synapses
# ------------------------------------------------------------------------------
cfg.synWeightFractionEE = [0.5, 0.5]  # E->E AMPA to NMDA ratio
cfg.synWeightFractionEI = [0.5, 0.5]  # E->I AMPA to NMDA ratio
cfg.synWeightFractionSOME = [0.9, 0.1]  # SOM -> E GABAASlow to GABAB ratio

cfg.synsperconn = {"HH_full": 5, "HH_reduced": 1, "HH_simple": 1, "PT": 3}
cfg.AMPATau2Factor = 1.0

# ------------------------------------------------------------------------------
# Network
# ------------------------------------------------------------------------------
cfg.singleCellPops = False  # Create pops with 1 single cell (to debug)
cfg.weightNorm = True  # use weight normalization
cfg.weightNormThreshold = 4.0  # weight normalization factor threshold

cfg.addConn = True
cfg.scale = 1.0
cfg.sizeY = 1350.0
cfg.sizeX = 300.0
cfg.sizeZ = 300.0
cfg.correctBorderThreshold = 150.0

cfg.L5BrecurrentFactor = 1.0
cfg.ITinterFactor = 1.0
cfg.strengthFactor = 1.0

cfg.EEGain = 0.5
cfg.EIGain = 1.0
cfg.IEGain = 1.0
cfg.IIGain = 1.0

cfg.IEdisynapticBias = None  # increase prob of I->Ey conns if Ex->I and Ex->Ey exist

# ------------------------------------------------------------------------------
## E->I gains
cfg.EPVGain = 1.0
cfg.ESOMGain = 1.0

# ------------------------------------------------------------------------------
## I->E gains
cfg.PVEGain = 1.0
cfg.SOMEGain = 1.0

# ------------------------------------------------------------------------------
## I->I gains
cfg.PVSOMGain = None  # 0.25
cfg.SOMPVGain = None  # 0.25
cfg.PVPVGain = None  # 0.75
cfg.SOMSOMGain = None  # 0.75

# ------------------------------------------------------------------------------
## I->E/I layer weights (L2/3+4, L5, L6)
cfg.IEweights = [0.8, 0.8, 1.0]
cfg.IIweights = [1.2, 1.0, 1.0]

cfg.IPTGain = 1.0
cfg.IFullGain = 1.0

# ------------------------------------------------------------------------------
# Subcellular distribution
# ------------------------------------------------------------------------------
cfg.addSubConn = True

# ------------------------------------------------------------------------------
# Long range inputs
# ------------------------------------------------------------------------------
cfg.addLongConn = True
cfg.numCellsLong = 1000  # num of cells per population
cfg.noiseLong = 1.0  # firing rate random noise
cfg.delayLong = 5.0  # (ms)
cfg.weightLong = 0.5  # corresponds to unitary connection somatic EPSP (mV)
cfg.startLong = 0  # start at 0 ms
# increaced motor thalamus TPO and TVL to 0-10Hz (from TPO [0, 5] and TVL [0, 2.5])
cfg.ratesLong = {
    "TPO": [0, 5],
    "TVL": [0, 2.5],
    "S1": [0, 5],
    "S2": [0, 5],
    "cM1": [0, 2.5],
    "M2": [0, 2.5],
    "OC": [0, 5],
}


## input pulses
cfg.addPulses = 0
# cfg.pulse = {'pop': 'None', 'start': 1000, 'end': 1200, 'rate': [0, 20], 'noise': 0.8}
# cfg.pulse2 = {'pop': 'None', 'start': 1000, 'end': 1200, 'rate': [0, 20], 'noise': 0.5, 'duration': 500}


# ------------------------------------------------------------------------------
# Current inputs
# ------------------------------------------------------------------------------
cfg.addIClamp = 0

cfg.IClamp1 = {
    "pop": "IT5B",
    "sec": "soma",
    "loc": 0.5,
    "start": 0,
    "dur": 1000,
    "amp": 0.50,
}


# ------------------------------------------------------------------------------
# NetStim inputs
# ------------------------------------------------------------------------------
cfg.addNetStim = False

## pop, sec, loc, synMech, start, interval, noise, number, weight, delay
cfg.NetStim1 = {
    "pop": "IT2",
    "ynorm": [0, 1],
    "sec": "soma",
    "loc": 0.5,
    "synMech": ["AMPA"],
    "synMechWeightFactor": [1.0],
    "start": 500,
    "interval": 1000.0 / 60.0,
    "noise": 0.0,
    "number": 60.0,
    "weight": 30.0,
    "delay": 0,
}
