import ROOT
from collections import defaultdict
from array import array
import numpy as np
import sys

file = ROOT.TFile("trackingNtuple.root")
tree = file.Get("tree")

# Your provided map
nmodules_map = {
    1: {  # VB
        0: [16, 5, 1],
        1: [16, 5, 1],
        2: [15, 5, 1],
        3: [15, 5, 1],
        4: [21, 5, 1],
        5: [21, 5, 1],
        6: [29, 5, 1],
        7: [29, 5, 1],
    },
    2: {  # VE
        0: [1, 16, 2],
        1: [1, 16, 2],
        2: [1, 16, 2],
        3: [1, 16, 2],
        4: [1, 16, 2],
        5: [1, 16, 2],
        6: [1, 16, 2],
        7: [1, 16, 2],
    },
}

observables = {
    "pt" : [0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 3.0, 5.0, 10, 15., 25, 50],
    "eta": np.linspace(-2.5, 2.5, 21).tolist(),
    "phi": np.linspace(-np.pi, np.pi, 21).tolist(),
    }

denom_cuts = {
    (1, 0) : {
        "pt" : {"pt": (0., 9999.), "abseta": (0., 0.8)},
        "eta": {"pt": (1., 9999.), "abseta": (0., 2.5)},
        "phi": {"pt": (1., 9999.), "abseta": (0., 0.8)},
        },
    (1, 1) : {
        "pt" : {"pt": (0., 9999.), "abseta": (0., 0.8)},
        "eta": {"pt": (1., 9999.), "abseta": (0., 2.5)},
        "phi": {"pt": (1., 9999.), "abseta": (0., 0.8)},
        },
    (1, 2) : {
        "pt" : {"pt": (0., 9999.), "abseta": (0., 0.8)},
        "eta": {"pt": (1., 9999.), "abseta": (0., 2.5)},
        "phi": {"pt": (1., 9999.), "abseta": (0., 0.8)},
        },
    (1, 3) : {
        "pt" : {"pt": (0., 9999.), "abseta": (0., 0.8)},
        "eta": {"pt": (1., 9999.), "abseta": (0., 2.5)},
        "phi": {"pt": (1., 9999.), "abseta": (0., 0.8)},
        },
    (1, 4) : {
        "pt" : {"pt": (0., 9999.), "abseta": (0., 0.5)},
        "eta": {"pt": (1., 9999.), "abseta": (0., 2.5)},
        "phi": {"pt": (1., 9999.), "abseta": (0., 0.5)},
        },
    (1, 5) : {
        "pt" : {"pt": (0., 9999.), "abseta": (0., 0.5)},
        "eta": {"pt": (1., 9999.), "abseta": (0., 2.5)},
        "phi": {"pt": (1., 9999.), "abseta": (0., 0.5)},
        },
    (1, 6) : {
        "pt" : {"pt": (0., 9999.), "abseta": (0., 0.5)},
        "eta": {"pt": (1., 9999.), "abseta": (0., 2.5)},
        "phi": {"pt": (1., 9999.), "abseta": (0., 0.5)},
        },
    (1, 7) : {
        "pt" : {"pt": (0., 9999.), "abseta": (0., 0.5)},
        "eta": {"pt": (1., 9999.), "abseta": (0., 2.5)},
        "phi": {"pt": (1., 9999.), "abseta": (0., 0.5)},
        },
    (2, 0) : {
        "pt" : {"pt": (0., 9999.), "abseta": (0., 1.8)},
        "eta": {"pt": (0., 9999.), "abseta": (0., 2.5)},
        "phi": {"pt": (0., 9999.), "abseta": (0., 1.8)},
        },
    (2, 1) : {
        "pt" : {"pt": (0., 9999.), "abseta": (0., 1.8)},
        "eta": {"pt": (0., 9999.), "abseta": (0., 2.5)},
        "phi": {"pt": (0., 9999.), "abseta": (0., 1.8)},
        },
    (2, 2) : {
        "pt" : {"pt": (0., 9999.), "abseta": (0., 1.8)},
        "eta": {"pt": (0., 9999.), "abseta": (0., 2.5)},
        "phi": {"pt": (0., 9999.), "abseta": (0., 1.8)},
        },
    (2, 3) : {
        "pt" : {"pt": (0., 9999.), "abseta": (0., 1.8)},
        "eta": {"pt": (0., 9999.), "abseta": (0., 2.5)},
        "phi": {"pt": (0., 9999.), "abseta": (0., 1.8)},
        },
    (2, 4) : {
        "pt" : {"pt": (0., 9999.), "abseta": (0., 1.8)},
        "eta": {"pt": (0., 9999.), "abseta": (0., 2.5)},
        "phi": {"pt": (0., 9999.), "abseta": (0., 1.8)},
        },
    (2, 5) : {
        "pt" : {"pt": (0., 9999.), "abseta": (0., 1.8)},
        "eta": {"pt": (0., 9999.), "abseta": (0., 2.5)},
        "phi": {"pt": (0., 9999.), "abseta": (0., 1.8)},
        },
    (2, 6) : {
        "pt" : {"pt": (0., 9999.), "abseta": (0., 1.8)},
        "eta": {"pt": (0., 9999.), "abseta": (0., 2.5)},
        "phi": {"pt": (0., 9999.), "abseta": (0., 1.8)},
        },
    (2, 7) : {
        "pt" : {"pt": (0., 9999.), "abseta": (0., 1.8)},
        "eta": {"pt": (0., 9999.), "abseta": (0., 2.5)},
        "phi": {"pt": (0., 9999.), "abseta": (0., 1.8)},
        },
    }

# Histogram containers
h_den = {}
h_num = {}

# Pre-book histograms only for valid system/layer combinations
for system in nmodules_map:
    for layer in nmodules_map[system]:
        h_den[(system, layer)] = {}
        h_num[(system, layer)] = {}
        for obs in observables:
            den = ROOT.TH1F(f"h_den_{obs}_sys{system}_layer{layer}", "", len(observables[obs])-1, array('f', observables[obs]))
            num = ROOT.TH1F(f"h_num_{obs}_sys{system}_layer{layer}", "", len(observables[obs])-1, array('f', observables[obs]))
            key = (system, layer)
            h_den[key][obs] = den
            h_num[key][obs] = num

# Loop over tree
for entry in tree:
    sim_ids = list(entry.sim_id)
    sim_pts = list(entry.sim_pt)
    sim_etas = list(entry.sim_eta)
    sim_phis = list(entry.sim_phi)
    sim_is_denom = list(entry.sim_isDenom)
    
    md_sim_ids = list(entry.md_sim_id)
    md_layers = list(entry.md_layer)
    md_systems = list(entry.md_system)

    # Map sim_id -> set of (system, layer) where it was reconstructed
    simid_to_syslayers = defaultdict(set)
    for i in range(len(md_sim_ids)):
        sid = int(md_sim_ids[i])
        system = int(md_systems[i])
        layer = int(md_layers[i])
        if system in nmodules_map and layer in nmodules_map[system]:
            simid_to_syslayers[sid].add((system, layer))

    # Loop over sim particles
    for i, sid in enumerate(sim_ids):
        sid = int(sid)
        pt = sim_pts[i]
        eta = sim_etas[i]
        abseta = abs(eta)
        phi = sim_phis[i]
        v = {"pt": pt, "eta": eta, "phi": phi, "abseta": abseta}

        # Fill histograms
        for system in nmodules_map:
            for layer in nmodules_map[system]:
                key = (system, layer)
                for obs in observables:
                    is_pass = True
                    for cutps in denom_cuts[key][obs]:
                        is_pass = is_pass and v[cutps] >= denom_cuts[key][obs][cutps][0]
                        is_pass = is_pass and v[cutps] <  denom_cuts[key][obs][cutps][1]
                    if is_pass:
                        h_den[key][obs].Fill(v[obs])
                        if key in simid_to_syslayers.get(sid, []):
                            h_num[key][obs].Fill(v[obs])
                        else:
                            print(f"key: {key} , v: {v} , ")

# Optional: save to file
outfile = ROOT.TFile("eff.root", "RECREATE")
for system in nmodules_map:
    for layer in nmodules_map[system]:
        key = (system, layer)
        for obs in observables:
            h_num[key][obs].Write()
            h_den[key][obs].Write()
outfile.Close()

