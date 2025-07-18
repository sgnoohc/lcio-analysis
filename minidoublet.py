from array import array
import os
from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1D, TH2D, TFile, TLorentzVector, TTree, TMath, TVector3, std, gInterpreter
from math import *
from optparse import OptionParser

# Enable ROOT's automatic C++ STL vector handling
gInterpreter.Declare("#include <vector>")

# VERBOSE=True
VERBOSE=False
STOPEVENT=10

def vprint(*args, **kwargs):
    if VERBOSE:
        print(*args, **kwargs)

#########################
# parameters

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile trackingNtuple.root',
                  type=str, default='trackingNtuple.root')
(options, args) = parser.parse_args()

tree = TTree("tree", "tree")

## Create branches
col_i = {}
col_i["nmd"] = array('i', [0])
col_i["nhit"] = array('i', [0])
col_i["nsim"] = array('i', [0])
for i_n in col_i:
    tree.Branch(i_n, col_i[i_n], f"{i_n}/I")
col_vf = {}
col_vf["hit_x"] = std.vector('float')()
col_vf["hit_y"] = std.vector('float')()
col_vf["hit_z"] = std.vector('float')()
col_vf["hit_system"] = std.vector('float')()
col_vf["hit_layer"] = std.vector('float')()
col_vf["hit_side"] = std.vector('float')()
col_vf["hit_module"] = std.vector('float')()
col_vf["hit_sensor"] = std.vector('float')()
col_vf["md_x_lower"] = std.vector('float')()
col_vf["md_y_lower"] = std.vector('float')()
col_vf["md_z_lower"] = std.vector('float')()
col_vf["md_x_upper"] = std.vector('float')()
col_vf["md_y_upper"] = std.vector('float')()
col_vf["md_z_upper"] = std.vector('float')()
col_vf["md_system"] = std.vector('float')()
col_vf["md_layer"] = std.vector('float')()
col_vf["md_side"] = std.vector('float')()
col_vf["md_module"] = std.vector('float')()
col_vf["md_sensor"] = std.vector('float')()
col_vf["md_sim_id"] = std.vector('float')()
col_vf["sim_id"] = std.vector('float')()
col_vf["sim_pt"] = std.vector('float')()
col_vf["sim_eta"] = std.vector('float')()
col_vf["sim_phi"] = std.vector('float')()
col_vf["sim_pdgId"] = std.vector('float')()
col_vf["sim_q"] = std.vector('float')()
col_vf["sim_isDenom"] = std.vector('float')()
col_vf["sim_hasmd"] = std.vector('float')()
for vf_n in col_vf:
    tree.Branch(vf_n, col_vf[vf_n])

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

tracker_systems = [
    "VB",
    "VE",
    # "IB",
    # "IE",
    # "OB",
    # "OE",
]

tracker_systems_map = {
    1: "VB",
    2: "VE",
    3: "IB",
    4: "IE",
    5: "OB",
    6: "OE",
}

nmodules_map = {
    1: { # VB
        0: [16, 5, 1], # 16 modules, 5 sensors, 1 sides
        1: [16, 5, 1], # 16 modules, 5 sensors, 1 sides
        2: [15, 5, 1], # 15 modules, 5 sensors, 1 sides
        3: [15, 5, 1], # 15 modules, 5 sensors, 1 sides
        4: [21, 5, 1], # 21 modules, 5 sensors, 1 sides
        5: [21, 5, 1], # 21 modules, 5 sensors, 1 sides
        6: [29, 5, 1], # 29 modules, 5 sensors, 1 sides
        7: [29, 5, 1], # 29 modules, 5 sensors, 1 sides
        },
    2: { # VE
        0: [1, 16, 2], # 0 modules, 16 sensors, 2 sides (-1 or 1)
        1: [1, 16, 2], # 0 modules, 16 sensors, 2 sides (-1 or 1)
        2: [1, 16, 2], # 0 modules, 16 sensors, 2 sides (-1 or 1)
        3: [1, 16, 2], # 0 modules, 16 sensors, 2 sides (-1 or 1)
        4: [1, 16, 2], # 0 modules, 16 sensors, 2 sides (-1 or 1)
        5: [1, 16, 2], # 0 modules, 16 sensors, 2 sides (-1 or 1)
        6: [1, 16, 2], # 0 modules, 16 sensors, 2 sides (-1 or 1)
        7: [1, 16, 2], # 0 modules, 16 sensors, 2 sides (-1 or 1)
        },
    }

module_keys = []
lower_module_keys = []
for system, layers in nmodules_map.items():
    for layer, (nmod, nsens, nsides) in layers.items():
        for module in range(nmod):
            for sensor in range(nsens):
                for side in [-1, 1] if nsides == 2 else [0]:
                    key = (system, layer, module, sensor, side)
                    module_keys.append(key)
                    if layer % 2 == 0:
                        lower_module_keys.append(key)

# Note: VX first layer
# VXBarrel
# Layer 0-1 = 0-15 modules, 0-4 sensors
# Layer 2-3 = 0-14 modules, 0-4 sensors
# Layer 4-5 = 0-20 modules, 0-4 sensors
# Layer 6-7 = 0-28 modules, 0-4 sensors
# VXEndcap
# Layer ALL = 0 modules, 0-15 sensors, -1 or 1 side

def dPhiThreshold(lhit_v3):
    Bfield = 5
    kRinv1GeVf = (2.99792458e-3 * Bfield)
    k2Rinv1GeVf = kRinv1GeVf / 2.
    sinAlphaMax = 0.95
    ptCut = 1
    rt = lhit_v3.Perp()
    miniSlope = asin(min(rt * k2Rinv1GeVf / ptCut, sinAlphaMax))
    error = 0
    return miniSlope + sqrt(error**2)

# loop over all events in the file
for ievent, event in enumerate(reader):

    if ievent == 0:
        for collection_name in event.getCollectionNames():
            col = event.getCollection(collection_name)
            coltype = str(col.getTypeName())
            print(f"{str(collection_name):50} {coltype:50}")
            if coltype == "LCRelation":
                relation = UTIL.LCRelationNavigator(col)
                print(f"{relation.getFromType()} -> {relation.getToType()}")

    if ievent % 100 == 0:
        print("Processing event " + str(ievent))

    # Clear all the branches
    for i_n in col_i:
        col_i[i_n][0] = -999
    for vf_n in col_vf:
        col_vf[vf_n].clear()

    hits = {}
    minidoublets = {}
    mcp_covered = []
    for module_key in module_keys:
        hits[module_key] = []
        minidoublets[module_key] = []

    hitCollections = {}
    hitRelationCollections = {}
    hitRelations = {}

    for tracker_system in tracker_systems:
        hitCollections[f"{tracker_system}"] = event.getCollection(f"{tracker_system}TrackerHits")
        hitRelationCollections[f"{tracker_system}"] = event.getCollection(f"{tracker_system}TrackerHitsRelations")
        hitRelations[f"{tracker_system}"] = UTIL.LCRelationNavigator(hitRelationCollections[f"{tracker_system}"])
        encoding = hitCollections[f"{tracker_system}"].getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder = UTIL.BitField64(encoding)
        for hit in hitCollections[f"{tracker_system}"]:
            cellID = int(hit.getCellID0())
            decoder.setValue(cellID)
            system = decoder["system"].value()
            layer = decoder["layer"].value()
            side = decoder["side"].value()
            module = decoder["module"].value()
            sensor = decoder["sensor"].value()
            hits[(system, layer, module, sensor, side)].append(hit)
            col_vf["hit_x"].push_back(hit.getPosition()[0])
            col_vf["hit_y"].push_back(hit.getPosition()[1])
            col_vf["hit_z"].push_back(hit.getPosition()[2])
            col_vf["hit_system"].push_back(system)
            col_vf["hit_layer"].push_back(layer)
            col_vf["hit_side"].push_back(side)
            col_vf["hit_module"].push_back(module)
            col_vf["hit_sensor"].push_back(sensor)

    if VERBOSE:
        print("hits")
        for key in hits:
            if len(hits[key]) > 0:
                print(f"key: {key} , hits[key]: {hits[key]} , ")

    for lower_module_key in lower_module_keys:

        upper_module_keys = []
        upper_module_keys.append((
            lower_module_key[0],
            lower_module_key[1] + 1,
            lower_module_key[2],
            lower_module_key[3],
            lower_module_key[4],
        ))

        if lower_module_key[0] == 1: # if VXB
            if lower_module_key[3] == 2: # if sensor = 2 then go both -z and +z
                upper_module_keys.append((
                    lower_module_key[0],
                    lower_module_key[1] + 1,
                    lower_module_key[2],
                    lower_module_key[3] - 1,
                    lower_module_key[4],
                ))
                upper_module_keys.append((
                    lower_module_key[0],
                    lower_module_key[1] + 1,
                    lower_module_key[2],
                    lower_module_key[3] + 1,
                    lower_module_key[4],
                ))
            elif lower_module_key[3] == 1: # if sensor = 1 then go -z
                upper_module_keys.append((
                    lower_module_key[0],
                    lower_module_key[1] + 1,
                    lower_module_key[2],
                    lower_module_key[3] - 1,
                    lower_module_key[4],
                ))
            elif lower_module_key[3] == 3: # if sensor = 3 then go +z
                upper_module_keys.append((
                    lower_module_key[0],
                    lower_module_key[1] + 1,
                    lower_module_key[2],
                    lower_module_key[3] + 1,
                    lower_module_key[4],
                ))

        for lhit in hits[lower_module_key]:
            for upper_module_key in upper_module_keys:
                for uhit in hits[upper_module_key]:
                    lhit_x = lhit.getPosition()[0]
                    lhit_y = lhit.getPosition()[1]
                    lhit_z = lhit.getPosition()[2]
                    uhit_x = uhit.getPosition()[0]
                    uhit_y = uhit.getPosition()[1]
                    uhit_z = uhit.getPosition()[2]

                    cellID = int(lhit.getCellID0())
                    decoder.setValue(cellID)
                    system = decoder["system"].value()
                    layer = decoder["layer"].value()
                    side = decoder["side"].value()
                    module = decoder["module"].value()
                    sensor = decoder["sensor"].value()

                    lhit_v3 = TVector3(lhit_x, lhit_y, lhit_z)
                    uhit_v3 = TVector3(uhit_x, uhit_y, uhit_z)

                    md_dz = lhit_v3.z() - uhit_v3.z()
                    md_dphi = lhit_v3.DeltaPhi(uhit_v3)
                    md_dphiChange = lhit_v3.DeltaPhi(uhit_v3 - lhit_v3)
                    md_dphiChange_threshold = dPhiThreshold(lhit_v3)

                    vprint(f"lhit_x: {lhit_x} , lhit_y: {lhit_y} , lhit_z: {lhit_z} , uhit_x: {uhit_x} , uhit_y: {uhit_y} , uhit_z: {uhit_z} , ")
                    vprint(f"md_dz: {md_dz} , md_dphi: {md_dphi} , md_dphiChange: {md_dphiChange} , md_dphiChange_threshold: {md_dphiChange_threshold} , ")

                    # dz cut
                    if abs(md_dz) > 5: continue

                    # absolute dphi cut
                    if abs(md_dphi) > md_dphiChange_threshold: continue

                    # dphi change cut
                    if abs(md_dphiChange) > md_dphiChange_threshold: continue

                    # Accept minidoublet
                    minidoublets[lower_module_key].append((lhit, uhit))

                    col_vf["md_x_lower"].push_back(lhit_x)
                    col_vf["md_y_lower"].push_back(lhit_y)
                    col_vf["md_z_lower"].push_back(lhit_z)
                    col_vf["md_x_upper"].push_back(uhit_x)
                    col_vf["md_y_upper"].push_back(uhit_y)
                    col_vf["md_z_upper"].push_back(uhit_z)
                    col_vf["md_system"].push_back(system)
                    col_vf["md_layer"].push_back(layer)
                    col_vf["md_side"].push_back(side)
                    col_vf["md_module"].push_back(module)
                    col_vf["md_sensor"].push_back(sensor)

                    # compute mcp_id
                    l_simhits = hitRelations[f"{tracker_systems_map[system]}"].getRelatedToObjects(lhit)
                    u_simhits = hitRelations[f"{tracker_systems_map[system]}"].getRelatedToObjects(uhit)
                    if len(l_simhits) > 0 and len(u_simhits) > 0:
                        l_mcp_id = l_simhits[0].getMCParticle().id()
                        u_mcp_id = u_simhits[0].getMCParticle().id()
                        if l_mcp_id == u_mcp_id:
                            mcp_covered.append(l_mcp_id)
                            col_vf["md_sim_id"].push_back(l_mcp_id)
                        else:
                            col_vf["md_sim_id"].push_back(-999)
                    else:
                        col_vf["md_sim_id"].push_back(-999)

    # set the counter
    col_i["nmd"][0] = col_vf["md_x_lower"].size()
    col_i["nhit"][0] = col_vf["hit_x"].size()

    # vprint(minidoublets)

    # get all MC particle that I want to track
    mcps = event.getCollection("MCParticle")
    for mcp in mcps:
        mcp_id = mcp.id()
        mcp_px = mcp.getMomentum()[0]
        mcp_py = mcp.getMomentum()[1]
        mcp_pz = mcp.getMomentum()[2]
        mcp_energy = mcp.getEnergy()
        tlv = TLorentzVector()
        tlv.SetPxPyPzE(mcp_px, mcp_py, mcp_pz, mcp_energy)
        pt = tlv.Pt()
        eta = tlv.Eta()
        phi = tlv.Phi()
        q = mcp.getCharge()
        pdgid = mcp.getPDG()
        isdenom = pt > 1 and abs(eta) < 1.8
        hasmd = mcp_id in mcp_covered
        vprint(f"mcp_id={mcp_id} pt={pt:.3f} eta={eta:.3f} phi={phi:.3f} pdgid={pdgid}")

        col_vf["sim_id"].push_back(mcp_id)
        col_vf["sim_pt"].push_back(pt)
        col_vf["sim_eta"].push_back(eta)
        col_vf["sim_phi"].push_back(phi)
        col_vf["sim_pdgId"].push_back(pdgid)
        col_vf["sim_q"].push_back(q)
        col_vf["sim_isDenom"].push_back(isdenom)
        col_vf["sim_hasmd"].push_back(hasmd)

    tree.Fill()

    if VERBOSE and ievent == STOPEVENT:
        break

reader.close()
output_file = TFile(options.outFile, 'RECREATE')
tree.Write()
output_file.Close()

# for tracker_subdet in tracker_subdets:
#     hitCollection = event.getCollection(f"{tracker_subdet}TrackerHits")
#     hitRelationCollection = event.getCollection(f"{tracker_subdet}TrackerHitsRelations")
#     hitRelation = UTIL.LCRelationNavigator(hitRelationCollection)
#     print(f"{hitRelation.getFromType()} -> {hitRelation.getToType()}")
#     hit_type_subdet = hit_type_map[tracker_subdet]

#     for hit in hitCollection:
#         hit_x[0] = hit.getPosition()[0]
#         hit_y[0] = hit.getPosition()[1]
#         hit_z[0] = hit.getPosition()[2]
#         hit_type[0] = hit_type_subdet
#         print(f"{hit.id()} {hit_x[0]:20} {hit_y[0]:20} {hit_z[0]:20} reco hit {hit.getU()[0]} {hit.getU()[1]} {hit.getdU()} {hit.getType()}")
#         simhits = hitRelation.getRelatedToObjects(hit)
#         nsimhit[0] = len(simhits)
#         for simhit in simhits:
#             simhit_x[0] = simhit.getPosition()[0]
#             simhit_y[0] = simhit.getPosition()[1]
#             simhit_z[0] = simhit.getPosition()[2]
#             print(f"{simhit.id()} {simhit_x[0]:20} {simhit_y[0]:20} {simhit_z[0]:20} matched simhit")
#             break
#         tree.Fill()

# simhitCollection = event.getCollection("VertexBarrelCollection")
# for simhit in simhitCollection:
#     simhit_x[0] = simhit.getPosition()[0]
#     simhit_y[0] = simhit.getPosition()[1]
#     simhit_z[0] = simhit.getPosition()[2]
#     mcp = simhit.getMCParticle()
#     mcp_px = mcp.getMomentum()[0]
#     mcp_py = mcp.getMomentum()[1]
#     mcp_pz = mcp.getMomentum()[2]
#     mcp_mass = mcp.getMass()
#     print(f"{simhit.id()} {simhit_x[0]:20} {simhit_y[0]:20} {simhit_z[0]:20} {simhit.getCellID0()} {simhit.getCellID1()} all simhit {mcp.getPDG()} {mcp_px} {mcp_py} {mcp_pz} {mcp_mass}")

# simhitCollection = event.getCollection("VertexBarrelCollectionConed")
# for simhit in simhitCollection:
#     simhit_x[0] = simhit.getPosition()[0]
#     simhit_y[0] = simhit.getPosition()[1]
#     simhit_z[0] = simhit.getPosition()[2]
#     print(f"{simhit.id()} {simhit_x[0]:20} {simhit_y[0]:20} {simhit_z[0]:20} all simhit")

# hitCol = event.getCollection("VBTrackerHitsConed")
# for hit in hitCol:
#     hit_x[0] = hit.getPosition()[0]
#     hit_y[0] = hit.getPosition()[1]
#     hit_z[0] = hit.getPosition()[2]
#     hit_type[0] = hit_type_subdet
#     print(f"{hit.id()} {hit_x[0]:20} {hit_y[0]:20} {hit_z[0]:20} reco hit ")





















# float SDL::CPU::MiniDoublet::dPhiThreshold(const SDL::CPU::Hit& lowerHit, const SDL::CPU::Module& module,const float dPhi, const float dz)
# {
#     // =================================================================
#     // Various constants
#     // =================================================================
#     const float kRinv1GeVf = (2.99792458e-3 * 3.8);
#     const float k2Rinv1GeVf = kRinv1GeVf / 2.;
#     // const float ptCut = PTCUT;
#     // const float sinAlphaMax = 0.95;
#     float ptCut = 1;
#     // std::cout <<  " module.layer(): " << module.layer() <<  std::endl;
#     // if (module.layer() == 6 or module.layer() == 5)
#     // {
#     //     ptCut = 0.96;
#     // }
#     float sinAlphaMax = 0.95;
#     // if (module.layer() == 6)
#     // {
#     //     sinAlphaMax = 2.95;
#     // }
#     // p2Sim.directionT-r2Sim.directionT smearing around the mean computed with ptSim,rSim
#     // (1 sigma based on 95.45% = 2sigma at 2 GeV)
#     std::array<float, 6> miniMulsPtScaleBarrel {0.0052, 0.0038, 0.0034, 0.0034, 0.0032, 0.0034};
#     std::array<float, 5> miniMulsPtScaleEndcap {0.006, 0.006, 0.006, 0.006, 0.006}; //inter/extra-polated from L11 and L13 both roughly 0.006 [larger R have smaller value by ~50%]
#     //mean of the horizontal layer position in y; treat this as R below
#     std::array<float, 6> miniRminMeanBarrel {21.8, 34.6, 49.6, 67.4, 87.6, 106.8}; // TODO: Update this with newest geometry
#     std::array<float, 5> miniRminMeanEndcap {131.4, 156.2, 185.6, 220.3, 261.5};// use z for endcaps // TODO: Update this with newest geometry

#     // =================================================================
#     // Computing some components that make up the cut threshold
#     // =================================================================
#     float rt = lowerHit.rt();
#     unsigned int iL = module.layer() - 1;
#     const float miniSlope = std::asin(std::min(rt * k2Rinv1GeVf / ptCut, sinAlphaMax));
#     const float rLayNominal = ((module.subdet() == SDL::CPU::Module::Barrel) ? miniRminMeanBarrel[iL] : miniRminMeanEndcap[iL]);
#     const float miniPVoff = 0.1 / rLayNominal;
#     const float miniMuls = ((module.subdet() == SDL::CPU::Module::Barrel) ? miniMulsPtScaleBarrel[iL] * 3.f / ptCut : miniMulsPtScaleEndcap[iL] * 3.f / ptCut);
#     const bool isTilted = module.subdet() == SDL::CPU::Module::Barrel and module.side() != SDL::CPU::Module::Center;
#     const bool tiltedOT123 = true;
#     const float pixelPSZpitch = 0.15;
#     const unsigned int detid = ((module.moduleLayerType() == SDL::CPU::Module::Pixel) ?  module.partnerDetId() : module.detId());
#     const float drdz = tiltedGeometry.getDrDz(detid);
#     const float miniTilt = ((isTilted && tiltedOT123) ? 0.5f * pixelPSZpitch * drdz / sqrt(1.f + drdz * drdz) / moduleGapSize(module) : 0);

#     // Compute luminous region requirement for endcap
#     const float deltaZLum = 15.f;
#     // const float miniLum = abs(dPhi * deltaZLum/dz); // Balaji's new error
#     const float miniLum = fabs(dPhi * deltaZLum/dz); // Balaji's new error
#     // const float miniLum = abs(deltaZLum / lowerHit.z()); // Old error


#     // =================================================================
#     // Return the threshold value
#     // =================================================================
#     // Following condition is met if the module is central and flatly lying
#     if (module.subdet() == SDL::CPU::Module::Barrel and module.side() == SDL::CPU::Module::Center)
#     {
#         return miniSlope + sqrt(pow(miniMuls, 2) + pow(miniPVoff, 2));
#     }
#     // Following condition is met if the module is central and tilted
#     else if (module.subdet() == SDL::CPU::Module::Barrel and module.side() != SDL::CPU::Module::Center) //all types of tilted modules
#     {
#         return miniSlope + sqrt(pow(miniMuls, 2) + pow(miniPVoff, 2) + pow(miniTilt * miniSlope, 2));
#     }
#     // If not barrel, it is Endcap
#     else
#     {
#         return miniSlope + sqrt(pow(miniMuls, 2) + pow(miniPVoff, 2) + pow(miniLum, 2));
#     }

# }

    # for tracker_system in tracker_systems:
    #     hitCollection = event.getCollection(f"{tracker_system}TrackerHits")
    #     hitRelationCollection = event.getCollection(f"{tracker_system}TrackerHitsRelations")
    #     hitRelation = UTIL.LCRelationNavigator(hitRelationCollection)
    #     encoding = hitCollection.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
    #     decoder = UTIL.BitField64(encoding)
    #     # <constant name="GlobalTrackerReadoutID"     type="string" value="system:5,side:-2,layer:6,module:11,sensor:8"/>
    #     # <constant name="GlobalCalorimeterReadoutID" type="string" value="system:5,side:-2,module:8,stave:4,layer:9,submodule:4,x:32:-16,y:-16"/>
    #     for hit in hitCollection:
    #         cellID = int(hit.getCellID0())
    #         decoder.setValue(cellID)
    #         _hit_layer[0] = decoder['layer'].value()
    #         _hit_system[0] = decoder["system"].value()
    #         _hit_side[0] = decoder["side"].value()
    #         _hit_module[0] = decoder["module"].value()
    #         _hit_sensor[0] = decoder["sensor"].value()
    #         _hit_x[0] = hit.getPosition()[0]
    #         _hit_y[0] = hit.getPosition()[1]
    #         _hit_z[0] = hit.getPosition()[2]
    #         simhits = hitRelation.getRelatedToObjects(hit)
    #         _simhit_x[0] = simhits[0].getPosition()[0]
    #         _simhit_y[0] = simhits[0].getPosition()[1]
    #         _simhit_z[0] = simhits[0].getPosition()[2]
    #         mcp = simhits[0].getMCParticle()
    #         mcp_px = mcp.getMomentum()[0]
    #         mcp_py = mcp.getMomentum()[1]
    #         mcp_pz = mcp.getMomentum()[2]
    #         mcp_energy = mcp.getEnergy()
    #         tlv = TLorentzVector()
    #         tlv.SetPxPyPzE(mcp_px, mcp_py, mcp_pz, mcp_energy)
    #         pt = tlv.Pt()
    #         eta = tlv.Eta()
    #         phi = tlv.Phi()
    #         pdgid = mcp.getPDG()
    #         mcp_id = mcp.id()

    #         vprint(
    #             f"tracker_system={tracker_system} "
    #             f"layer={_hit_layer[0]} "
    #             f"system={_hit_system[0]} "
    #             f"side={_hit_side[0]} "
    #             f"module={_hit_module[0]} "
    #             f"sensor={_hit_sensor[0]} "
    #             f"hit_x={_hit_x[0]:.3f} "
    #             f"hit_y={_hit_y[0]:.3f} "
    #             f"hit_z={_hit_z[0]:.3f}, "
    #             f"pt={pt:.3f}, "
    #             f"eta={eta:.3f}, "
    #             f"phi={phi:.3f} "
    #             f"pdgid={pdgid} "
    #             f"mcp_id={mcp_id}"
    #         )
    #         tree.Fill()
